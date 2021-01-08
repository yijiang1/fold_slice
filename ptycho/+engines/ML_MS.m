% [ p, fdb ] = ML_MS( p )

% Academic License Agreement
%
% Source Code
%
% Introduction 
% •	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
%   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the cSAXS 
%   ptychography MATLAB package computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").
%
% Terms and Conditions of the LICENSE
% 1.	LICENSOR grants to LICENSEE a royalty-free, non-exclusive license to use the PROGRAM for academic, non-commercial purposes, upon the terms and conditions 
%       hereinafter set out and until termination of this license as set forth below.
% 2.	LICENSEE acknowledges that the PROGRAM is a research tool still in the development stage. The PROGRAM is provided without any related services, improvements 
%       or warranties from LICENSOR and that the LICENSE is entered into in order to enable others to utilize the PROGRAM in their academic activities. It is the 
%       LICENSEE’s responsibility to ensure its proper use and the correctness of the results.”
% 3.	THE PROGRAM IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR 
%       A PARTICULAR PURPOSE AND NONINFRINGEMENT OF ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. IN NO EVENT SHALL THE LICENSOR, THE AUTHORS OR THE COPYRIGHT 
%       HOLDERS BE LIABLE FOR ANY CLAIM, DIRECT, INDIRECT OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION WITH THE PROGRAM OR THE USE 
%       OF THE PROGRAM OR OTHER DEALINGS IN THE PROGRAM.
% 4.	LICENSEE agrees that it will use the PROGRAM and any modifications, improvements, or derivatives of PROGRAM that LICENSEE may create (collectively, 
%       "IMPROVEMENTS") solely for academic, non-commercial purposes and that any copy of PROGRAM or derivatives thereof shall be distributed only under the same 
%       license as PROGRAM. The terms "academic, non-commercial", as used in this Agreement, mean academic or other scholarly research which (a) is not undertaken for 
%       profit, or (b) is not intended to produce works, services, or data for commercial use, or (c) is neither conducted, nor funded, by a person or an entity engaged 
%       in the commercial use, application or exploitation of works similar to the PROGRAM.
% 5.	LICENSEE agrees that it shall make the following acknowledgement in any publication resulting from the use of the PROGRAM or any translation of the code into 
%       another computing language:
%       "Data processing was carried out using the cSAXS ptychography MATLAB package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
%       Scherrer Institut, Switzerland."
%
% Additionally, any publication using the package, or any translation of the code into another computing language should cite for difference map:
% P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379–382 (2008). 
%   (doi: 10.1126/science.1158573),
% for maximum likelihood:
% P. Thibault and M. Guizar-Sicairos, Maximum-likelihood refinement for coherent diffractive imaging, New J. Phys. 14, 063004 (2012). 
%   (doi: 10.1088/1367-2630/14/6/063004),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
%   (doi: 10.1364/OE.24.029089).
% 6.	Except for the above-mentioned acknowledgment, LICENSEE shall not use the PROGRAM title or the names or logos of LICENSOR, nor any adaptation thereof, nor the 
%       names of any of its employees or laboratories, in any advertising, promotional or sales material without prior written consent obtained from LICENSOR in each case.
% 7.	Ownership of all rights, including copyright in the PROGRAM and in any material associated therewith, shall at all times remain with LICENSOR, and LICENSEE 
%       agrees to preserve same. LICENSEE agrees not to use any portion of the PROGRAM or of any IMPROVEMENTS in any machine-readable form outside the PROGRAM, nor to 
%       make any copies except for its internal use, without prior written consent of LICENSOR. LICENSEE agrees to place the following copyright notice on any such copies: 
%       © All rights reserved. PAUL SCHERRER INSTITUT, Switzerland, Laboratory for Macromolecules and Bioimaging, 2017. 
% 8.	The LICENSE shall not be construed to confer any rights upon LICENSEE by implication or otherwise except as specifically set forth herein.
% 9.	DISCLAIMER: LICENSEE shall be aware that Phase Focus Limited of Sheffield, UK has an international portfolio of patents and pending applications which relate 
%       to ptychography and that the PROGRAM may be capable of being used in circumstances which may fall within the claims of one or more of the Phase Focus patents, 
%       in particular of patent with international application number PCT/GB2005/001464. The LICENSOR explicitly declares not to indemnify the users of the software 
%       in case Phase Focus or any other third party will open a legal action against the LICENSEE due to the use of the program.
% 10.	This Agreement shall be governed by the material laws of Switzerland and any dispute arising out of this Agreement or use of the PROGRAM shall be brought before 
%       the courts of Zürich, Switzerland. 

function [ p, fdb ] = ML_MS( p )
import utils.*
    global opt_time
    fdb.status = core.engine_status;
    
    opt_time = 0;
    verbose(1, 'Starting multi-slice non-linear optimization')
    
    % ===== 3ML =====
    N_layer = p.N_layer;
    Ny = p.asize(1);
    Nx = p.asize(2);
    lambda = p.lambda; 
    k  = 2*pi/lambda; 
    
    % ----- Calculate the propagator
    [Xp,Yp] = meshgrid(([1:p.asize(2)]-floor(p.asize(2)/2)+1)*p.dx_spec(2), ([1:p.asize(1)]-floor(p.asize(1)/2)+1)*p.dx_spec(1));
    Xp = ifftshift(Xp); 
    Yp = ifftshift(Yp);
    dx = Xp(1,2)-Xp(1,1);
    Fx = Xp/(Nx*dx^2);
    dy = Yp(2,1)-Yp(1,1);
    Fy = Yp/(Ny*dy^2);
%     for n = 1:N_layer-1
%         propagation{n} = exp( 1j*k*p.delta_z(n)*sqrt( 1-(lambda*Fx).^2-(lambda*Fy).^2 ) );
%         propagation_back{n} = exp( 1j*k*(-p.delta_z(n))*sqrt( 1-(lambda*Fx).^2-(lambda*Fy).^2 ) );
%     end
    p.Fx = Fx;
    p.Fy = Fy;
    
    % ----- Initialization      
    if ~isfield(p,'object_layers') && ndims(p.object{1}) < 4
        if (length(p.ms_init_ob_fraction) ~= p. N_layer) || sum(p.ms_init_ob_fraction)~=1
            verbose(0,'-- (Initialization) p.ms_init_ob_fraction not good, will use 1/N_layer for all layers');
            p.ms_init_ob_fraction = ones(1,p. N_layer)/p. N_layer;
        end
    
        for obnum = 1:p.numobjs
            if (isfield(p,'initial_iterate_object') && strcmp(p.initial_iterate_object,'file')) || (max(angle(p.object{obnum}(:)))-min(angle(p.object{obnum}(:)))) > 1.5*pi % specify input or propagating results from another engine
                ob_phase{obnum} = engines.ML_MS.fun_ramp_unwrap(p.object{obnum}, p.asize);
            else
                ob_phase{obnum} = angle(p.object{obnum});
            end
            for n = 1:N_layer            
                for obmode = 1:p.object_modes            
                    object_layer{obnum}{obmode}{n} = abs(p.object{obnum}).^(p.ms_init_ob_fraction(n)) .* exp(1i*ob_phase{obnum}.*p.ms_init_ob_fraction(n));
                end

                if p.use_display
                    figure(100); 
                    subplot(N_layer,2,2*n-1); imagesc(abs(object_layer{1}{1}{n})); colormap bone; axis equal xy tight; caxis([0 2]); colorbar; drawnow;
                    subplot(N_layer,2,2*n); imagesc(angle(object_layer{1}{1}{n})); colormap bone; axis equal xy tight; caxis([-pi pi]); colorbar; drawnow;
                    if n==1
                        title('Initial image, layer 1');
                    end
                end
            end
        end	
    elseif ndims(p.object{1}) == 4
        verbose(0,'-- (Initialization) Using previous multilayer results');
        for obnum = 1:p.numobjs
            for n = 1:N_layer  
                for obmode = 1:p.object_modes            
                    object_layer{obnum}{obmode}{n} = double(p.object{obnum}(:,:,obmode,n)); 
                end
            end
        end
        
    else
        verbose(0,'-- (Initialization) Using previous MS results');
        for n = 1:N_layer
           for obnum = 1:p.numobjs   
               for obmode = 1:p.object_modes   
                    object_layer{obnum}{obmode}{n} = p.object_layers{n};
               end
           end
        end
    end
    probes = p.probes;
        
    recon_time_tic = tic;
    recon_time = [];
    delta_z_iter = [];   
    % ========== 
    
    if p.probe_mask_bool
        
        if p.probe_mask_use_auto
            verbose(2, 'Using a probe mask from probe autocorrelation.');
            to_threshold = -real(auto);
        else
            verbose(2, 'Using a circular probe mask.');
            [x,y] = meshgrid(-p.asize(2)/2:floor((p.asize(2)-1)/2),-p.asize(1)/2:floor((p.asize(1)-1)/2));
            to_threshold = (x.^2 + y.^2);
            clear x y
        end
        to_threshold_flat = reshape(to_threshold, [prod(p.asize) 1]);
        [~, ind] = sort(to_threshold_flat);
        probe_mask_flat = zeros([prod(p.asize) 1]);
        probe_mask_flat(ind(1:ceil(p.probe_mask_area * prod(p.asize)))) = 1;
        p.probe_mask = reshape(probe_mask_flat, p.asize);
        clear to_threshold to_threshold_flat dummy ind probe_mask_flat
    else
        p.probe_mask = ones(p.asize);
    end
    
    % Taking care to pass some needed functions in p
    fnorm = sqrt(prod(p.asize));

    %%% Optimization error metric
    if isfield(p,'opt_errmetric'),
        switch lower(p.opt_errmetric)
            case 'l1'
                verbose(1, 'Using ML-L1 error metric'),
            case 'l2'
                verbose(1,'Using ML-L2 error metric'),
            case 'poisson'
                verbose(1,'Using ML-Poisson'),
            otherwise
                error([p.opt_errmetric ' is not defined'])
                return;
        end
    else
        p.opt_errmetric = 'poisson';
        verbose(1, 'Using default Poisson error metric')
    end
    
    %%% Set specific variables needed for different metrics %%%
    switch lower(p.opt_errmetric)
        case 'poisson'
            fmag2 = p.fmag.^2;
            fmag2renorm = fmag2/p.renorm^2;
            initialerror = p.renorm^2*sum( p.fmask(:).*( (fmag2renorm(:)+0.5).*log(fmag2renorm(:)+1) ...
                - fmag2renorm(:) ...
                - 1 + 0.5*log(2*pi) + 1./(12*(fmag2renorm(:)+1))  ...
                - 1./(360*(fmag2renorm(:)+1).^3) + 1./(1260*(fmag2renorm(:)+1).^5) )) ... %% Approximation to log(n!) http://www.johndcook.com/blog/2010/08/16/how-to-compute-log-factorial/
                + sum( fmag2(:)*log(renorm^2) );
            clear fmag2renorm
        case 'l1'
            initialerror = 0;
            fmag2 = 0;
        case 'l2'
            initialerror = 0;
            fmag2 = p.fmag.^2;
        otherwise
            error(['Error metric ' p.opt_errmetric 'is not defined'])
    end
    
    
    %%% Regularization
    Npix = 0;
    if p. reg_mu > 0
        for obnum = 1:p.numobjs
            Npix = Npix + p.object_size(obnum,1)*p.object_size(obnum,2);
        end
        Nm = prod(p.asize)*size(p.fmag,3);
        K = 8*Npix^2/(Nm*p.Nphot);
        creg = p.renorm^2*p.reg_mu/K;
    else
        creg = 0;
    end
    
    %%% Sieves preconditioning
    if any(p.smooth_gradient) ~= 0
        if length(p.smooth_gradient) <= 1
            % Hanning regularization
            auxi = fract_hanning_pad(512,512,0);
            auxi = fftshift(ifft2(auxi));
            smooth_gradient = real(auxi(256:258,256:258));  % Regularization kernel ( = 0 to omit)
        end
    else
        smooth_gradient = 0;
    end
    
    %% ===== 3ML main =====
    p.ms_opt_flags_local = p.ms_opt_flags;  
    
    if p.ms_opt_flags(3)
        N_iter_outer = ceil(p.ms_opt_iter/p.ms_opt_z_param(1)) + floor(p.ms_opt_iter/p.ms_opt_z_param(1));
    else
        N_iter_outer = 1;
        N_iter_inner = p.ms_opt_iter;
    end
    
    core.errorplot; % clear the persistent variable
    for iter_outer = 1:N_iter_outer
        
        if p.ms_opt_flags(3)         
            p.ms_opt_flags_local(3) = ~mod(iter_outer,2); % Alternates between updating delta_z, [0 1 0 1...]
            if p.ms_opt_flags_local(3)
                N_iter_inner = p. ms_opt_z_param(2); % when update delta_z
            else
                N_iter_inner = p. ms_opt_z_param(1); 
            end
        end
       
        optimize_object_layer = p.ms_opt_flags_local(1);
        optimize_probes = p.ms_opt_flags_local(2);
        optimize_delta_z = p.ms_opt_flags_local(3);
        delta_z_iter = [delta_z_iter; p.delta_z(:)'];
        
        % -- Arranging optimization vector 
        xopt = [];
        if optimize_object_layer
            for obnum = 1:p.numobjs
                for obmode = 1:p.object_modes
                    for n = 1:N_layer
                        xopt = [xopt; reshape([real(object_layer{obnum}{obmode}{n}(:)).'; imag(object_layer{obnum}{obmode}{n}(:)).'], [], 1)];
                    end
                end
            end
        else
            p.object_layer = object_layer;
        end
        if optimize_probes
            xopt = [xopt; reshape([real(probes(:)).'; imag(probes(:)).'], [], 1)];
        else
            p.probes = probes;
        end
        if optimize_delta_z
            xopt = [xopt; p.delta_z];
        end   
      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Main optimization loop %%%
        opt_time = tic;
        [tmp, p] = engines.ML.cgmin1('engines.ML_MS.gradient_ptycho_MS', xopt, N_iter_inner, p.opt_ftol, p.opt_xtol,...
            p, p.fmag, fmag2, p.fmask, p.numobjs, p.object_size, p.numprobs,...
            p.numscans, p.scanindexrange, initialerror, fnorm, p.probe_mask, p.plot.errtitlestring, [],...
            [], creg, smooth_gradient);

        opt_time = toc(opt_time);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
   
        % --- Error and time
        ms_opt_error = core.errorplot([]);   % Only reads the persistent variable
        ms_error(2, iter_outer) = ms_opt_error(end);
        recon_time(iter_outer) = toc(recon_time_tic);
        
        fprintf('[iter_outer %d, with flags %d %d %d] opt_time = %.2f min, total %.0f min\n\n', ...
            iter_outer, p.ms_opt_flags_local(1), p.ms_opt_flags_local(2), p.ms_opt_flags_local(3), opt_time/60, recon_time(end)/60);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Arrange solution vector %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if optimize_object_layer
            for obnum = 1:p.numobjs
                for obmode = 1:p.object_modes
                    for n = 1:N_layer
                        o_size = p.object_size(obnum, :);
                        o_numel = prod(o_size);
                        object_layer{obnum}{obmode}{n} = reshape(tmp(1:2:2*o_numel), o_size) + ...
                            1i*reshape(tmp(2:2:2*o_numel), o_size);
                        tmp = tmp(2*o_numel+1:end);
                    end
                end
            end
        end
        if optimize_probes
            probeelements = [p.asize p.numprobs p.probe_modes];
            probes = reshape(tmp(1:2:2*prod(probeelements)), probeelements) + ...
                1i*reshape(tmp(2:2:2*prod(probeelements)), probeelements);
            tmp = tmp(2*prod(probeelements)+1:end);
        end    
        if optimize_delta_z
            delta_z = tmp;
            tmp = tmp(size(delta_z)+1:end);
            fprintf(' Optimized delta_z = %.4f um\n',delta_z(1)*1e6);
            if delta_z<0
                fprintf(' Force delta_z = 0 um\n');
                p.delta_z = 0;
            else
                p.delta_z = delta_z;
            end
        end  
        
        if ~isempty(tmp)
            warning('Temporary vector is not empty, optimized values not assigned');
        end
    
    end % end iter_outer
        
    verbose(2, 'Finished');
    verbose(2, 'Time elapsed in optimization refinement: %f seconds', opt_time);
    
    p.delta_z_iter = delta_z_iter;
    p.error_metric.iteration = 1:length(ms_opt_error); 
    p.error_metric.value = ms_opt_error; 
    p.error_metric.err_metric = p.opt_errmetric; 
    p.error_metric.method = p.name;
    p.recon_time = recon_time;
    
    for obnum = 1:p.numobjs
        for n = 1:N_layer  
            for obmode = 1:p.object_modes            
                p.object{obnum}(:,:,obmode,n) = object_layer{obnum}{obmode}{n}; 
            end
        end
    end
        
               
    %%%%%%%%%%%%%%%%%
    %%% Last plot %%%
    %%%%%%%%%%%%%%%%%
    if p.use_display||p.store_images
        p.plot.extratitlestring = sprintf(' (%dx%d) - ML', p.asize(2), p.asize(1));
        core.analysis.plot_results(p, 'use_display', p.use_display, 'store_images', p.store_images);
    end
    core.errorplot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% end optimization refinement %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



