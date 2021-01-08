%ML Maximum Likelihood refinement
%
%   Publications most relevant to the Maximum Likelihood refinement
%       + M. Guizar-Sicairos and J. R. Fienup, "Phase retrieval with transverse
%       translation diversity: a nonlinear optimization approach," Opt. Express 16, 7264-7278 (2008)
%       + P. Thibault and M. Guizar-Sicairos, "Maximum-likelihood refinement for
%       coherent diffractive imaging," New J. Phys. 14, 063004 (2012).

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

function [ p, fdb ] = ML( p )

    import utils.*

    global opt_time
    fdb.status = core.engine_status;
    core.errorplot;

    opt_time = 0;
    verbose(3, 'Starting non-linear optimization')
    
    if p.probe_mask_bool
        
        if p.probe_mask_use_auto
            verbose(3, 'Using a probe mask from probe autocorrelation.');
            to_threshold = -real(auto);
        else
            verbose(3, 'Using a circular probe mask.');
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

    % Arranging optimization vector
    xopt = []; % Optimization vector
    if p.opt_flags(1) == 1
        for obnum = 1:p.numobjs
            xopt = [xopt; real(p.object{obnum}(:)); imag(p.object{obnum}(:))];
        end
    end
    if p.opt_flags(2) == 1,
        xopt = [xopt; real(p.probes(:)); imag(p.probes(:))];
    end
    
    xopt = single(xopt);  % assumed by the MEX scripts
    p.fmag = single(p.fmag); 
    
    
    % if flags(3) == 1,
    %     xopt = [xopt;x;y];
    % else
    %     fixed.x = x;
    %     fixed.y = y;
    % end
    if isempty(xopt),
        error('At least one element of flags must be 1'),
    end
    
    %%% Optimization error metric
    if isfield(p,'opt_errmetric'),
        switch lower(p.opt_errmetric)
            case 'l1'
                verbose(2, 'Using ML-L1 error metric'),
            case 'l2'
                verbose(2,'Using ML-L2 error metric'),
            case 'poisson'
                verbose(2,'Using ML-Poisson'),
            otherwise
                error([p.opt_errmetric ' is not defined'])
                return;
        end
    else
        p.opt_errmetric = 'poisson';
        verbose(2, 'Using default Poisson error metric')
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
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Main optimization loop %%%
    opt_time = tic;
    [tmp,p] = engines.ML.cgmin1('engines.ML.gradient_ptycho',xopt,p.opt_iter,p.opt_ftol,p.opt_xtol,p,fmag2,initialerror,fnorm,...
        creg, smooth_gradient); % ob and probes are passed to use as fixed variables when p.opt_flags is zero
 
    opt_time = toc(opt_time);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Arrange solution vector %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if p.opt_flags(1) == 1,
        for obnum = 1:p.numobjs
            objectelements = [p.object_size(obnum,1) p.object_size(obnum,2) p.object_modes];
            p.object{obnum} = reshape(tmp(1:prod(objectelements)), objectelements) + ...
                1i*reshape(tmp(prod(objectelements)+1:2*prod(objectelements)), objectelements);
            tmp = tmp(2*prod(objectelements)+1:end);
        end
    end
    if p.opt_flags(2) == 1,
        probeelements = [p.asize p.numprobs p.probe_modes];
        p.probes = reshape(tmp(1:prod(probeelements)),probeelements) + ...
            1i*reshape(tmp(prod(probeelements)+1:2*prod(probeelements)),probeelements);
        tmp = tmp(2*prod(probeelements)+1:end);
    end
    %     if flags(3) == 1,
    %         x = tmp(1:params.n);
    %         y = tmp(params.n+1:2*params.n);
    %     end
    if ~isempty(tmp)
        warning('Temporary vector is not empty, optimized values not assigned');
    end
    
    verbose(3, 'Finished');
    verbose(3, 'Time elapsed in optimization refinement: %f seconds', opt_time);
    
    %%%%%%%%%%%%%%%%%
    %%% Last plot %%%
    %%%%%%%%%%%%%%%%%
    if p.use_display||p.save.store_images
        p.plot.extratitlestring = sprintf(' (%dx%d) - ML', p.asize(2), p.asize(1));
        core.analysis.plot_results(p, 'use_display', p.use_display);
    end
    core.errorplot;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% end optimization refinement %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



