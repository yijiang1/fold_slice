% [ p, fdb ] = DM_MS( p )

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

function [ p, fdb ] = DM_MS( p )
import utils.verbose

global proj1_time  objproj_time probeproj_time elsewheretime proj2_time plot_time

fdb.status = core.engine_status;

% ===== 3DM  =====
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
for n = 1:N_layer-1
    propagation{n} = exp( 1j*k*p.delta_z(n)*sqrt( 1-(lambda*Fx).^2-(lambda*Fy).^2 ) );
    propagation_back{n} = exp( 1j*k*(-p.delta_z(n))*sqrt( 1-(lambda*Fx).^2-(lambda*Fy).^2 ) );
end
p.Fx = Fx;
p.Fy = Fy;
    
% ----- Initialization
if (length(p.ms_init_ob_fraction) ~= p. N_layer) || sum(p.ms_init_ob_fraction)~=1
    fprintf('-- (Initialization) p.ms_init_ob_fraction bad, will use 1/N_layer for all layers \n');
    p.ms_init_ob_fraction = ones(1,p. N_layer)/p. N_layer;
end

for obnum = 1:p.numobjs
    if (isfield(p,'initial_iterate_object') && strcmp(p.initial_iterate_object,'file')) || (max(angle(p.object{obnum}(:)))-min(angle(p.object{obnum}(:)))) > 1.5*pi % specify input or propagating results from another engine
        ob_phase{obnum} = engines.ML_MS.fun_ramp_unwrap(p.object{obnum}, p.asize);
    else
        ob_phase{obnum} = angle(p.object{obnum});
    end
    for n = 1:N_layer            
        for obmode = 1:p.object_modes % object modes?                
            object_layer{obnum}{obmode}{n} = abs(p.object{obnum}).^(p.ms_init_ob_fraction(n)) .* exp(1i*ob_phase{obnum}.*p.ms_init_ob_fraction(n));
        end

        figure(100); 
        subplot(N_layer,2,2*n-1); imagesc(abs(object_layer{1}{1}{n})); colormap bone; axis equal xy tight; caxis([0 2]); colorbar; drawnow;
        subplot(N_layer,2,2*n); imagesc(angle(object_layer{1}{1}{n})); colormap bone; axis equal xy tight; caxis([-pi pi]); colorbar; drawnow;
        if n==1
            title('Initial image, layer 1');
        end
    end
end	
probes = p.probes;

recon_time_tic = tic;
recon_time = [];
delta_z_iter = [];   
% ========== 

% mex or matlab - object_update / probe_update / Fourier_loop
if ~isfield(p, 'use_mex')
    p.use_mex = zeros(1,3);
elseif size(p.use_mex,1) == 1
    if any(p.use_mex)
        verbose(3, 'Using mex files for DM.')
    end
    p.use_mex = repmat(p.use_mex,1,3);
end


% Starting iterate
dsize = [p.asize, p.numpos];
% iter is a large array and is only needed for Matlab Difference map
%iter = cell(p.numscans,1);
% fmag = cell(p.numscans,1);
% for ii = 1:p.numscans
%     %iter{ii} = complex(zeros([dsize p.object_modes*p.probe_modes]));
%     fmag{ii} = double(p.fmag(:,:,p.scanidxs{ii})); % can cause memory duplication if it was not double, migrate to single in future !!!
% end
if p.object_modes == 1 && (p.numscans == 1 ||  length(p.scanidxs{1}) == length(p.scanidxs{1}) )
    obj_proj = complex(zeros([p.asize, length(p.scanidxs{1})]));
end


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
        

% check what is the actual size of the object 
for i = 1:p.numscans
    p.object_size(i,:) = size(p.object{i});
    ob{i} = double(p.object{i});
end
p.probes = double(p.probes);

for obnum = 1:p.numobjs
    avob{obnum} = zeros([p.object_size(obnum,:) p.object_modes]);
end

% get views from probes and objects
% for ii = 1:p.numscans
%     for prmode = 1:p.probe_modes
%         for obmode = 1:p.object_modes
%             iter_mode_ind = prmode+(obmode-1)*p.probe_modes;
%             iter(:,:,p.scanidxs{ii},iter_mode_ind) = bsxfun(@times, p.probes(:,:,p.share_probe_ID(ii),prmode), core.get_projections(p, p.object{p.share_object_ID(ii)}, ii));
%         end
%     end
% end 
% for ii = 1:p.numscans
%     if p.share_object
%         obnum = 1;
%     else
%         obnum = ii;
%     end
%     if p.object_modes == 1
%         % faster version without extra memory allocation 
%         obj_proj = core.get_projections(p, ob{obnum}, ii,obj_proj);
%         iter{ii} = bsxfun(@times, p.probes, obj_proj);
%     else
%         for obmode = 1:p.object_modes
%             obj_proj = core.get_projections(p, ob{obnum}(:,:,obmode), ii);
%             iter_mode_ind = (1:p.probe_modes)+(obmode-1)*p.probe_modes;
%             iter{ii}(:,:,:,iter_mode_ind) = bsxfun(@times, p.probes, obj_proj);
%         end
%     end
% end


err = nan(p.number_iterations,1);
rfact = nan(p.number_iterations,1);


% A power bound (relaxed Fourier) that scales with number of photons per diffraction pattern
p.power_bound = p.count_bound*p.renorm^2;

% Set indices for user supplied flat mask  p.object_flat_region
if ~isempty(p.object_flat_region)
    p.userflatregion = true;
    if (p.numscans>1)&&(~p.share_object)
        error('Object flat region not yet implemented for multiple objects. Set object_flat_region = []')
    end
    if any(p.object_size ~= size(p.object_flat_region))
        error('Mask p.object_flat_region does not match size of object');
    else
        p.userflatind = find(p.object_flat_region == 1);% Find indices of flat region
    end
else
    p.userflatregion = false;
end
    
for obnum = 1:p.numobjs
    avob{obnum} = zeros([p.object_size(obnum,:) p.object_modes]);
end

numav = 0;
cfact = p.probe_regularization *p.numpts;
if p.share_probe
    cfact = sum(cfact);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 3D Main Difference map loop %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   


% Prepare statistics
proj1_time = 0;
proj2_time = 0;
plot_time = 0;
objproj_time =   0;
probeproj_time = 0;
elsewheretime =  0;


cfact = p.probe_regularization *p.numpts;
if p.share_probe
    cfact = sum(cfact);
end

numav = 0;

%%%% Started off from the copy from multi-slice branch
positions = p.positions;
scanindexrange = p.scanindexrange;
whichtokeep = [1:(size(p.fmask,3))];

a2 = p.asize(1) * p.asize(2);
obmode = 1;
prmode = 1;
ms_error = core.errorplot; % Clears the persistent variable
    
for it = 1:p.number_iterations

    % ------ Fourier projection of view_layer{N_layer} ------
    fprintf(' ---- Fourier modulus constraint ---- \n');

    p1 = zeros(p.asize)+eps*(1+1i);
    p2 = zeros(p.asize)+eps*(1+1i);
    f  = zeros(p.asize)+eps*(1+1i);
    ph = zeros(p.asize)+eps*(1+1i);
    df = zeros(p.asize)+eps*(1+1i);
    af = zeros(p.asize);
    fdev = zeros(p.asize);
    fdev2 = zeros(p.asize);
    fmaski = zeros(p.asize);
    rf = 0;
    rf_nrm = 0;

    er2 = 0;
    func_L1 = 0;
    
    for ii = 1:p.numscans        
        
        if p.share_object
            obnum = 1;
        else
            obnum = ii;
        end
        if p.share_probe
            prnum = 1;
        else
            prnum = ii;
        end
        fmask_per_scan = ndims(p.fmask) == 3; 

        if it==1
            fprintf(' Calculate view from initial guess.. \n');
            probe_layer{prnum}{prmode}{1} = probes(:,:,prnum,prmode);

            for jj=scanindexrange(ii,1):scanindexrange(ii,2)
                Indy = round(positions(jj,1)) + (1:p.asize(1));
                Indx = round(positions(jj,2)) + (1:p.asize(2));
                for n = 1:N_layer
                    if n==1
                        probe_use = probe_layer{prnum}{prmode}{1};
                    else
                        probe_use = probe_layer{prnum}{prmode}{n}{jj};
                    end
                    view_layer{n}{jj}  = probe_use .* object_layer{obnum}{obmode}{n}(Indy, Indx);
                    if n+1 <= N_layer
                        probe_layer{prnum}{prmode}{n+1}{jj} = ifft2(fft2(view_layer{n}{jj}) .* propagation{n});
                    end
                end
            end
            probe_layer_init = probe_layer;
            object_layer_init = object_layer;
            view_layer_init = view_layer;
        end

        fnorm = sqrt(a2);
        rf = 0;
        rf2 = 0;
        rf_nrm = 0;
        if ~fmask_per_scan
            fmaski = p.fmask;
        end

        fprintf(' Fourier projection.. scan = %d\n', ii);
        for jj=scanindexrange(ii,1):scanindexrange(ii,2)
            if fmask_per_scan
                fmaski = p.fmask(:,:,jj);
                fmaski = repmat(fmaski,[1 1 p.probe_modes*p.object_modes]);
            end
            Indy = round(positions(jj,1)) + (1:p.asize(1));
            Indx = round(positions(jj,2)) + (1:p.asize(2));
            Iq = 0;
            for obmode = 1:p.object_modes
                for prmode = 1:p.probe_modes
                    if N_layer==1
                        probe_use = probe_layer{prnum}{prmode}{1};
                    else
                        if ~isempty( find(~(whichtokeep-jj),1) ) % in the update ROI
                            probe_use = probe_layer{prnum}{prmode}{N_layer}{jj};
                        else
                            probe_use = probe_layer_init{prnum}{prmode}{N_layer}{jj};
                        end
                    end
                    iter_mode_ind = prmode+(obmode-1) * p.probe_modes;
                    p1(:,:,iter_mode_ind) = probe_use .* object_layer{obnum}{obmode}{N_layer}(Indy, Indx);
                    psiq = fft2(p1(:,:,iter_mode_ind))/fnorm;
                    Iq = Iq + abs(psiq).^2;
                end
            end
            f = fft2( 2*p1 - squeeze(view_layer{N_layer}{jj}) )/fnorm;
            af = abs(f); % Amplitude of f before projection
            ph = f ./ (af+1e-10);
            % Target of a perfect projection
            %   fmag_target = repmat(fmag(:,:,jj),[1 1
            %   p.probe_modes*p.object_modes]); % Naive target
            fmag_target = af.*repmat(p.fmag(:,:,jj)./sqrt(sum(af.^2,3)),[1 1 p.probe_modes*p.object_modes]);
            fdev = af - fmag_target;
            fdev2 = fmaski.*fdev.^2;
            power = sum(fdev2(:))/a2;
            if power > p.power_bound
                renorm = sqrt(p.power_bound / power);
                af = af.*(1-fmaski) + fmaski.*(fmag_target + fdev * renorm);
            end
            p2 = fnorm*ifft2(af .* ph);

            df = p2 - p1;
            view_layer{N_layer}{jj} = view_layer{N_layer}{jj} + reshape(df, size(view_layer{N_layer}{jj}));

            er2 = er2 + sum(abs(df(:)).^2);
            if p.compute_rfact
                rf = rf + sum(sum(abs( abs(fft2(p1)/fnorm) - p.fmag(:,:,jj) )));
                rf2 = rf2 + sum(sum( ( abs(fft2(p1)/fnorm) - p.fmag(:,:,jj) ).^2));
                rf_nrm = rf_nrm + sum(sum(p.fmag(:,:,jj)));
            end

            % Calculate L1 error 
            Fq = sqrt(Iq);
            if p.inv_intensity 
                alpha = sum(sum(fmaski.*p.fmag(:,:,jj).*Fq))/sum(sum(fmaski.*Fq.^2));
            else
                alpha = 1;
            end
            func_L1 = func_L1 + sum(sum( fmaski.*( alpha*Fq - p.fmag(:,:,jj)).^2 ));
        end
    end

    % ------ Object and probe updates 
    fprintf(' ---- Object and probe updates ---- \n');
    for n = N_layer:(-1):1
        for obnum = 1:p.numobjs
            object_layer_old{obnum}{obmode}{n} = object_layer{obnum}{obmode}{n};
            object_denom{obnum} = 1e-8 * ones(p.object_size(obnum,:));
            object_layer{obnum}{obmode}{n} = 1e-8*(1+1i)*ones([p.object_size(obnum,:)  p.object_modes]);
        end

        fprintf(' Object update.. n = %d\n', n);
        for ii = 1:p.numscans
            if p.share_object
                obnum = 1;
            else
                obnum = ii;
            end
            if p.share_probe
                prnum = 1;
            else
                prnum = ii;
            end

            objtic = tic;
            for jj=scanindexrange(ii,1):scanindexrange(ii,2)
                Indy = round(positions(jj,1)) + (1:p.asize(1));
                Indx = round(positions(jj,2)) + (1:p.asize(2));
                for prmode = 1:p.probe_modes
                    if n==1
                        probe_use = probe_layer{prnum}{prmode}{1};
                    else
                        if ~isempty( find(~(whichtokeep-jj),1) ) % in the update ROI
                            probe_use = probe_layer{prnum}{prmode}{N_layer}{jj};
                        else
                            probe_use = probe_layer_init{prnum}{prmode}{N_layer}{jj};
                        end
                    end
                    if ~isempty( find(~(whichtokeep-jj),1) ) % in the update ROI
                        view_use = view_layer{n}{jj};
                    else
                        view_use = view_layer_init{n}{jj};
                    end
                    for obmode = 1:p.object_modes
                        iter_mode_ind = prmode+(obmode-1)*p.probe_modes;
                        object_layer{obnum}{obmode}{n}(Indy,Indx) = object_layer{obnum}{obmode}{n}(Indy,Indx) + conj(probe_use) .* view_use;
                    end
                    object_denom{obnum}(Indy,Indx) = object_denom{obnum}(Indy,Indx) + abs(probe_use).^2;
                end
            end
            objproj_time = objproj_time + toc(objtic);
        end % end ii

        for obnum = 1:p.numobjs
            object_layer{obnum}{obmode}{n} = object_layer{obnum}{obmode}{n} ./ object_denom{obnum};
            elsewheretic = tic;
            if p.clip_object
                aob = abs(object_layer{obnum}{obmode}{n});
                temp_phase = object_layer{obnum}{obmode}{n} ./ aob;
                indx_max = aob > p.clip_max;
                indx_min = aob < p.clip_min;
                object_layer{obnum}{obmode}{n}(indx_max) = p.clip_max*temp_phase(indx_max);
                object_layer{obnum}{obmode}{n}(indx_min) = p.clip_min*temp_phase(indx_min);
            end
            elsewheretime = elsewheretime + toc(elsewheretic);
        end

        fprintf(' Probe update.. n = %d\n', n);
        if n==1
            probe_layer_old{prnum}{prmode}{1} = probe_layer{prnum}{prmode}{1};
            for prnum = 1:p.numprobs
                probe_layer{prnum}{prmode}{1} = probe_layer{prnum}{prmode}{1} *cfact(prnum);
                probe_demon(:,:,prnum) = ones(p.asize)*cfact(prnum);
            end
        end
        for ii = 1:p.numscans
            if p.share_object
                obnum = 1;
            else
                obnum = ii;
            end
            if p.share_probe
                prnum = 1;
            else
                prnum = ii;
            end
            probetic = tic;

            for jj=scanindexrange(ii,1):scanindexrange(ii,2)
                Indy = round(positions(jj,1)) + (1:p.asize(1));
                Indx = round(positions(jj,2)) + (1:p.asize(2));
                if ~isempty( find(~(whichtokeep-jj),1) ) % in the update ROI
                    for obmode = 1:p.object_modes
                        for prmode = 1:p.probe_modes
                            iter_mode_ind = prmode+(obmode-1)*p.probe_modes;
                            if n>1
                                probe_layer{prnum}{prmode}{n}{jj} = conj(object_layer{obnum}{obmode}{n}(Indy,Indx)) .* view_layer{n}{jj} ...
                                    ./ abs(object_layer{obnum}{obmode}{n}(Indy,Indx)).^2;
                            else
                                probe_layer{prnum}{prmode}{1} = probe_layer{prnum}{prmode}{1} + ...
                                    conj(object_layer{obnum}{obmode}{1}(Indy,Indx)) .* view_layer{1}{jj};
                            end
                        end
                        % pr_denom = pr_denom + abs(ob{obnum}(Indy,Indx,obmode)).^2;
                        if n==1
                            probe_demon(:,:,prnum) = probe_demon(:,:,prnum) + abs(object_layer{obnum}{obmode}{1}(Indy,Indx)).^2;
                        end
                    end
                end
            end
        end

        if it >= p.probe_change_start
            if n==1
                for prnum = 1:p.numprobs
                    for prmode = 1:p.probe_modes
                        if p.probe_mask
                            probe_layer{prnum}{prmode}{1} = probe_mask .* probe_layer{prnum}{prmode}{1} ./ probe_demon(:,:,prnum);
                        else
                            probe_layer{prnum}{prmode}{1} = probe_layer{prnum}{prmode}{1} ./ probe_demon(:,:,prnum);
                        end
                    end
                end
            end
        else
            fprintf('           Probe0 NOT updated \n');
            probe_layer{prnum}{prmode}{1} = probes(:,:,prnum,prmode);
        end

        % --- Propagate reconstructed probe_layer to the previous view
        if n>=2
            for ii = 1:p.numscans
                if p.share_object
                    obnum = 1;
                else
                    obnum = ii;
                end
                if p.share_probe
                    prnum = 1;
                else
                    prnum = ii;
                end
                for jj=scanindexrange(ii,1):scanindexrange(ii,2)
                    if ~isempty( find(~(whichtokeep-jj),1) ) % in the update ROI
                        view_layer{n-1}{jj} = ifft2(fft2(probe_layer{prnum}{prmode}{n}{jj}) .* propagation_back{n-1});
                    end
                end
            end
        end

        figure(202); clf
        for nn = 1:N_layer
            subplot(N_layer,1,nn); 
            img0 = angle(object_layer_init{obnum}{obmode}{nn});
            img1 = angle(object_layer{obnum}{obmode}{nn});
            imagesc(img0-img1); axis xy equal tight; colorbar; colormap bone
            caxis([-0.5 0.5])
            title(sprintf('slice %d difference',nn));
        end
        drawnow
        %pause

    end % end n

    % ------ Temporarily stores object and probe (of the original order)
    for n = 1:N_layer
        for obnum = 1:p.numobjs
            object_layer_temp{obnum}{obmode}{n} = object_layer{obnum}{obmode}{n};
        end
    end
    for prnum = 1:p.numprobs
        probe_layer_temp{prnum}{prmode}{1} = probe_layer{prnum}{prmode}{1};
    end

    % ====== Reverse order: Update object
    if ~isempty(p.ms_reverse_order_iter) && (it >= p.ms_reverse_order_iter(1)) && (it <= p.ms_reverse_order_iter(2))
        fprintf(' Calculate view and probe (using objects and probe from the previous iteration).. \n');
        
        object_layer = object_layer_old;
        for ii = 1:p.numscans
            if p.share_object
                obnum = 1;
            else
                obnum = ii;
            end
            if p.share_probe
                prnum = 1;
            else
                prnum = ii;
            end
            probe_layer{prnum}{prmode}{1} = probe_layer_old{prnum}{prmode}{1};

%                 for jj=scanindexrange(ii,1):scanindexrange(ii,2)
%                     Indy = positions(jj,1) + (1:asize(1));
%                     Indx = positions(jj,2) + (1:asize(2));
%                     for n = 1:N_layer
%                         if n==1
%                             probe_use = probe_layer{prnum}{prmode}{1};
%                         else
%                             probe_use = probe_layer{prnum}{prmode}{n}{jj};
%                         end
%                         view_layer{n}{jj}  = probe_use .* object_layer{obnum}{obmode}{n}(Indy,Indx);
%                         if n+1 <= N_layer
%                             probe_layer{prnum}{prmode}{n+1}{jj} = ifft2(fft2(view_layer{n}{jj}) .* propagation{n});
%                         end
%                     end
%                 end
        end

        for n = 1:N_layer
            for obnum = 1:p.numobjs
                object_denom{obnum} = 1e-8 * ones(p.object_size(obnum,:));
                object_layer{obnum}{obmode}{n} = 1e-8*(1+1i)*ones([p.object_size(obnum,:)  p.object_modes]);
            end

            fprintf(' (Reverse Order) Object update.. n = %d\n', n);
            for ii = 1:p.numscans
                if p.share_object
                    obnum = 1;
                else
                    obnum = ii;
                end
                if p.share_probe
                    prnum = 1;
                else
                    prnum = ii;
                end

                objtic = tic;
                for jj=scanindexrange(ii,1):scanindexrange(ii,2)
                    Indy = positions(jj,1) + (1:p.asize(1));
                    Indx = positions(jj,2) + (1:p.asize(2));
                    for prmode = 1:p.probe_modes
                        if n==1
                            probe_use = probe_layer{prnum}{prmode}{1};
                        else
                            probe_use = probe_layer{prnum}{prmode}{n}{jj};
                        end
                        % --- Back propagate the measurement from N_layer
                        for layer = N_layer:(-1):n+1
                            probe_layer{prnum}{prmode}{layer}{jj} = conj(object_layer{obnum}{obmode}{layer}(Indy,Indx)) .* view_layer{layer}{jj} ...
                                    ./ abs(object_layer{obnum}{obmode}{layer}(Indy,Indx) + 1e-8).^2;
                            view_layer{layer-1}{jj} = ifft2(fft2(probe_layer{prnum}{prmode}{layer}{jj}) .* propagation_back{layer-1});
                        end
                        % --- Update object
                        for obmode = 1:p.object_modes
                            iter_mode_ind = prmode+(obmode-1)*p.probe_modes;
                            object_layer{obnum}{obmode}{n}(Indy,Indx) = object_layer{obnum}{obmode}{n}(Indy,Indx) + conj(probe_use) .* view_layer{n}{jj};
                        end
                        object_denom{obnum}(Indy,Indx) = object_denom{obnum}(Indy,Indx) + abs(probe_use).^2;
                    end
                end 
                objproj_time = objproj_time + toc(objtic);    
            end % end ii

            for obnum = 1:p.numobjs
                object_layer{obnum}{obmode}{n} = object_layer{obnum}{obmode}{n} ./ object_denom{obnum};
                elsewheretic = tic;
                if p.clip_object
                    aob = abs(object_layer{obnum}{obmode}{n});
                    temp_phase = object_layer{obnum}{obmode}{n} ./ aob;
                    indx_max = aob > p.clip_max;
                    indx_min = aob < p.clip_min;
                    object_layer{obnum}{obmode}{n}(indx_max) = p.clip_max*temp_phase(indx_max);
                    object_layer{obnum}{obmode}{n}(indx_min) = p.clip_min*temp_phase(indx_min);
                end
                elsewheretime = elsewheretime + toc(elsewheretic);

            end

            if n==1
                fprintf(' (Reverse Order) Probe update, only for n == 1\n');
                for prnum = 1:p.numprobs
                    probe_layer{prnum}{prmode}{1} = probe_layer{prnum}{prmode}{1} *cfact(prnum);
                    probe_demon(:,:,prnum) = ones(p.asize)*cfact(prnum);
                end

                for ii = 1:p.numscans
                    if p.share_object
                        obnum = 1;
                    else
                        obnum = ii;
                    end
                    if p.share_probe
                        prnum = 1;
                    else
                        prnum = ii;
                    end

                    probetic = tic;

                    for jj=scanindexrange(ii,1):scanindexrange(ii,2)
                        Indy = positions(jj,1) + (1:p.asize(1));
                        Indx = positions(jj,2) + (1:p.asize(2));
                        for obmode = 1:p.object_modes
                            for prmode = 1:p.probe_modes
                                iter_mode_ind = prmode+(obmode-1)*p.probe_modes;
                                probe_layer{prnum}{prmode}{1} = probe_layer{prnum}{prmode}{1} + ...
                                        conj(object_layer{obnum}{obmode}{1}(Indy,Indx)) .* view_layer{1}{jj};
                            end
                            % pr_denom = pr_denom + abs(ob{obnum}(Indy,Indx,obmode)).^2;
                            if n==1
                                probe_demon(:,:,prnum) = probe_demon(:,:,prnum) + abs(object_layer{obnum}{obmode}{1}(Indy,Indx)).^2;
                            end
                        end
                    end
                end

                for prnum = 1:p.numprobs
                    for prmode = 1:p.probe_modes
                        if p.probe_mask
                            probe_layer{prnum}{prmode}{1} = probe_mask .* probe_layer{prnum}{prmode}{1} ./ probe_demon(:,:,prnum);
                        else
                            probe_layer{prnum}{prmode}{1} = probe_layer{prnum}{prmode}{1} ./ probe_demon(:,:,prnum);
                        end
                    end
                end
            end

            % --- Calculate view{n} and probe{n+1} ----
            if n<=N_layer-1
                for ii = 1:p.numscans
                    if p.share_object
                        obnum = 1;
                    else
                        obnum = ii;
                    end
                    if p.share_probe
                        prnum = 1;
                    else
                        prnum = ii;
                    end
                    for jj=scanindexrange(ii,1):scanindexrange(ii,2)
                        Indy = positions(jj,1) + (1:p.asize(1));
                        Indx = positions(jj,2) + (1:p.asize(2));
                        if n==1
                            probe_use = probe_layer{prnum}{prmode}{1};
                        else
                            probe_use = probe_layer{prnum}{prmode}{n}{jj};
                        end
                        view_layer{n}{jj}  = probe_use .* object_layer{obnum}{obmode}{n}(Indy,Indx);
                        probe_layer{prnum}{prmode}{n+1}{jj} = ifft2(fft2(view_layer{n}{jj}) .* propagation{n});
                    end
                end
            end

        end % end n

        % ---- TAKE AVERAGE ----
        for n = 1:N_layer
            for obnum = 1:p.numobjs
                object_layer{obnum}{obmode}{n} = p.ratio_reverse*object_layer{obnum}{obmode}{n} + (1-p.ratio_reverse)*object_layer_temp{obnum}{obmode}{n};
            end
        end
        for prnum = 1:p.numprobs
            probe_layer{prnum}{prmode}{1} = p.ratio_reverse*probe_layer{prnum}{prmode}{1} + (1-p.ratio_reverse)*probe_layer_temp{prnum}{prmode}{1};
        end        

    end

    % ====== Forward model - Calculate views and probes for all layers
    fprintf(' ---- Forward model: Calculate views and probes for all layers.. \n');
    for ii = 1:p.numscans
        if p.share_object
            obnum = 1;
        else
            obnum = ii;
        end
        if p.share_probe
            prnum = 1;
        else
            prnum = ii;
        end
        for jj=scanindexrange(ii,1):scanindexrange(ii,2)
            Indy = round(positions(jj,1)) + (1:p.asize(1));
            Indx = round(positions(jj,2)) + (1:p.asize(2));
            for n = 1:N_layer
                if n==1
                    probe_use = probe_layer{prnum}{prmode}{1};
                else
                    if ~isempty( find(~(whichtokeep-jj),1) ) % in the update ROI
                        probe_use = probe_layer{prnum}{prmode}{N_layer}{jj};
                    else
                        probe_use = probe_layer_init{prnum}{prmode}{N_layer}{jj};
                    end
                end                    
                if ~isempty( find(~(whichtokeep-jj),1) ) % in the update ROI
                    view_layer{n}{jj}  = probe_use .* object_layer{obnum}{obmode}{n}(Indy,Indx);
                    if n+1 <= N_layer
                        probe_layer{prnum}{prmode}{n+1}{jj} = ifft2(fft2(view_layer{n}{jj}) .* propagation{n});
                    end
                end
            end
        end
    end

    if p.center_probe % Find probe centroid
        fprintf(' ##### p.center_probe: feature not implemented. #####\n');
    end


    % ------ Error metric
    ms_error(1,it) = sqrt(er2/(a2*sum(p.numpts))); % err(it) = sqrt(er2/(a2*sum(numpts)));
    ms_error(2,it) = func_L1;
    recon_time(it) = toc(recon_time_tic);
    delta_z_iter = [delta_z_iter; p.delta_z(:)'];
    if p.compute_rfact
        rfact(1,it) = rf/rf_nrm;
        rfact(2, it) = rf2/rf_nrm; % unused, can enter another error metric
        verbose(2, '[iter %d] Error: %.3f; R-factor: %.3f %%', it, ms_error(1,it), 100*rfact(1,it));
    else
        verbose(2,'[iter %d] Error: %.5f; func_L1: %.5f', it, ms_error(1,it), func_L1);
    end


   
end   

p.delta_z_iter = delta_z_iter;
p.error_metric.iteration = 1:size(ms_error,2); 
p.error_metric.value = ms_error(2,:);  
p.recon_time = recon_time;
p.error_metric.method = p.name;
p.error_metric.err_metric = p.opt_errmetric; 

for ii = 1:p.numscans
    if p.share_object
        obnum = 1;
    else
        obnum = ii;
    end
    if p.share_probe
        prnum = 1;
    else
        prnum = ii;
    end
    object = ones(size(object_layer{obnum}{obmode}{n})); 
    for n = 1:N_layer 
        object = object .* object_layer{obnum}{obmode}{n};   % Combine all layers
        p.object_layers{n} = object_layer{obnum}{obmode}{n}; % For each scan (object number) 
    end   
    p.object{obnum} = object;
    p.probes = probe_layer{prnum}{prmode}{1}; %probes(:,:,prnum,:);
end

%%%%%%%%%%%%%%%%%
%%% Last plot %%%
%%%%%%%%%%%%%%%%%
if p.use_display||p.store_images
    p.plot.extratitlestring = sprintf(' (%dx%d) - 3DM', p.asize(2), p.asize(1));
    core.analysis.plot_results(p, 'use_display', p.use_display, 'store_images', p.store_images);
end
core.errorplot;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% end main difference map loop %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


verbose(2, 'Finished difference map');
verbose(2, 'Time elapsed in projection 1: %f seconds', proj1_time);
verbose(2, '             in object projection: %f seconds', objproj_time);
verbose(2, '             in probe projection: %f seconds', probeproj_time);
%     verbose(2, '             in handpicked line: %f seconds', elsewheretime);
verbose(2, 'Time elapsed in projection 2: %f seconds', proj2_time);
verbose(2, 'Time spent plotting: %f seconds', plot_time);

% Average
for obnum = 1:p.numobjs
    if numav > 0
        avob{obnum} =  avob{obnum}/ numav;
    else
        avob{obnum} = ob{obnum};
    end
end
% R-factor
av_rfact = 0;
av_rfact_nrm = 0;
for ii = 1:p.numscans
    if p.share_object
        obnum = 1;
    else
        obnum = ii;
    end
    if p.share_probe
        prnum = 1;
    else
        prnum = ii;
    end
    for i=p.scanidxs{ii}
        Indy = round(p.positions(i,1)) + (1:p.asize(1));
        Indx = round(p.positions(i,2)) + (1:p.asize(2));
        fcalc = abs(fftn(avob{obnum}(Indy,Indx).*p.probes(:,:,prnum)))/sqrt(prod(p.asize));
        av_rfact = av_rfact + sum(sum(abs(fcalc - p.fmag(:,:,i))));
        av_rfact_nrm = av_rfact_nrm + sum(sum(p.fmag(:,:,i)));
    end
end
av_rfact = av_rfact / av_rfact_nrm;

% p.object = avob;


% save additional feedback into engine's fdb variable
fdb.rfact = rfact;
fdb.av_rfact = av_rfact;




end