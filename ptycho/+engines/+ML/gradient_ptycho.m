% Main code to compute error metric and gradient
% Jan 09 2013

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

function [func, grad, p] = gradient_ptycho(xopt,p,fmag2, initialerror,fnorm,creg,smooth_gradient)

import utils.verbose

%%% Initialize variables %%%
func = 0;    % Should be zero except for poisson (factorial factor)
for ii = 1:p.numobjs
    grado{ii} = zeros([p.object_size(ii,:) p.object_modes], 'like', xopt)+1i*eps;
end
gradp = zeros(p.asize(1),p.asize(2),p.numprobs,p.probe_modes, 'like', xopt)+1i*eps;
% gradx = zeros(n,1);
% grady = zeros(n,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Arrange optimization variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if p.opt_flags(1) == 1,
    for obnum = 1:p.numobjs
%         ob{obnum} = reshape(xopt(1:p.object_size(obnum,1)*p.object_size(obnum,2)),...
%             p.object_size(obnum,1),p.object_size(obnum,2)) + ...
%             1i*reshape(xopt(p.object_size(obnum,1)*p.object_size(obnum,2)+1:2*p.object_size(obnum,1)*p.object_size(obnum,2)),...
%             p.object_size(obnum,1),p.object_size(obnum,2));
        ob{obnum} = reshape(xopt(1:numel(grado{obnum})),...
            size(grado{obnum})) + ...
                 1i*reshape(xopt(numel(grado{obnum})+1:2*numel(grado{obnum})),...
            size(grado{obnum}));
        xopt = xopt(2*numel(grado{obnum})+1:end);
    end
end
if p.opt_flags(2) == 1,
    probes = reshape(xopt(1:numel(gradp)),size(gradp)) + ...
          1i*reshape(xopt(numel(gradp)+1:2*numel(gradp)),size(gradp));
    xopt = xopt(2*numel(gradp)+1:end);
end
%     if flags(3) == 1,
%         x = tmp(1:params.n);
%         y = tmp(params.n+1:2*params.n);
%     end

%%%%%%%%%%%%%%%%%%%%%
%%% Support error %%%
%%%%%%%%%%%%%%%%%%%%%
% Option to add later for a smooth support constraint error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute error metric %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Still to implement here: 
% + L1 and L2 metrics
% + Insensitive to multiplicative scale
if nargout > 1
    verbose(3,'Computing gradient')
end


for ii = 1:p.numscans
    prnum = p.share_probe_ID(ii); 
    obnum = p.share_object_ID(ii); 
    probe = probes(:,:,prnum,:); 
    Iq_all = 0;
    for obmode = 1:p.object_modes
        obj_proj{obmode} = core.get_projections(p, ob{obnum}(:,:,obmode), ii);
        psiq_all{obmode} = fft2(bsxfun(@times,obj_proj{obmode},probe/fnorm)); % view in Fourier domain
        Iq_all = Iq_all + sum(abs(psiq_all{obmode}).^2,4);
    end
        
    
    % Use implicit matlab paralelization to avoid computational overhead ,
    % currenly implemented only for L1 norm 
    if strcmpi(p.opt_errmetric,'l1')
        fmag = p.fmag(:,:,p.scanidxs{ii});
        fmask = p.fmask(:,:,p.scanidxs{ii});
        Fq = sqrt(Iq_all);
        %%% Invariant to intensity fluctuations
        if p.inv_intensity
            alpha = sum(sum(fmask.*fmag.*Fq))./sum(sum(fmag.*Fq.^2));
        else
            alpha = 1; 
        end
        func = sum(sum(sum(fmask.*( alpha.*Fq - fmag ).^2)));
        if nargout > 1 % Compute gradients
            for obmode = 1:p.object_modes
                chir = alpha.*ifft2(fmask.*( alpha - fmag./(Fq+eps) ).*psiq_all{obmode})*fnorm; % May not be needed for position optimization
                if p.opt_flags(1) == 1,
                    grado{obnum}(:,:,obmode) = core.set_projections(p, grado{obnum}(:,:,obmode), sum(2*conj(probe).*chir,4), ii);
                end
                if p.opt_flags(2) == 1
                    gradp(:,:,prnum,:) = gradp(:,:,prnum,:) ...
                        + sum(2*conj(obj_proj{obmode}).*chir,3);
                end
            end
        end
    else
      for jj = p.scanidxs{ii}  % Loop through diffraction patterns

            Indy = round(p.positions(jj,1)) + (1:p.asize(1));
            Indx = round(p.positions(jj,2)) + (1:p.asize(2));
            Iq = Iq_all(:,:,jj-p.scanidxs{ii}(1)+1);



            switch lower(p.opt_errmetric)
                    case 'poisson'
                        %%% Invariant to intensity fluctuations
                        if p.inv_intensity
                            % The numerator could be computed once outside
                            alpha = sum(sum(p.fmask(:,:,jj).*fmag2(:,:,jj)))/sum(sum(p.fmask(:,:,jj).*Iq));
                        else
                            alpha = 1;
                        end
                        func = func - sum(sum(p.fmask(:,:,jj).*( fmag2(:,:,jj).*log(alpha*Iq) - alpha*Iq )));

                        if nargout > 1 % Compute gradients
                            for obmode = 1:p.object_modes
                                psiq = psiq_all{obmode}(:,:,jj); % view in Fourier domain
                                chir = ifft2(p.fmask(:,:,jj).*( alpha - fmag2(:,:,jj)./Iq ).*psiq)*fnorm; % May not be needed for position optimization
                                for prmode = 1:p.probe_modes

                                    if p.opt_flags(1) == 1
                                        grado{obnum}(Indy,Indx,obmode) = grado{obnum}(Indy,Indx,obmode) ...
                                            + sum(2*conj(probe).*chir, 4);
                                    end
                                    if p.opt_flags(2) == 1
                                        gradp(:,:,prnum,:) = gradp(:,:,prnum,:) ...
                                            + 2*conj(ob{obnum}(Indy,Indx,obmode)).*chir;
                                    end                    
                                end
                            end
                        end
                case 'l2'
                    %%% Invariant to intensity fluctuations
                    if p.inv_intensity
                        alpha = sum(sum(p.fmask(:,:,jj).*fmag2(:,:,jj).*Iq))/sum(sum(p.fmask(:,:,jj).*Iq.^2));
                    else
                        alpha = 1;
                    end
                    tmp = alpha*Iq - fmag2(:,:,jj);
                    func = func + sum(sum(p.fmask(:,:,jj).*( tmp ).^2));
                    if nargout > 1 % Compute gradients
                        for obmode = 1:p.object_modes
                            psiq = psiq_all{obmode}(:,:,jj); % view in Fourier domain
                            chir = alpha*ifft2(2*p.fmask(:,:,jj).*( tmp ).*psiq)*fnorm; % May not be needed for position optimization
                            if p.opt_flags(1) == 1
                                grado{obnum}(Indy,Indx,obmode) = grado{obnum}(Indy,Indx,obmode) ...
                                    + sum(2*conj(probe).*chir,4);
                            end
                            if p.opt_flags(2) == 1
                                gradp(:,:,prnum,:) = gradp(:,:,prnum,:) ...
                                    + 2*conj(ob{obnum}(Indy,Indx,obmode)).*chir;
                            end
                        end
                    end

               case 'l1'
    %                 %%% Invariant to intensity fluctuations
    %                 if p.inv_intensity
    %                     alpha = sum(sum(p.fmask(:,:,jj).*p.fmag(:,:,jj).*Fq))/sum(sum(p.fmask(:,:,jj).*Fq.^2));
    %                 else
    %                     alpha = 1;
    %                 end
    %                 func = func + sum(sum(p.fmask(:,:,jj).*( alpha*Fq - p.fmag(:,:,jj) ).^2));
    %                 if nargout > 1 % Compute gradients
    %                     for obmode = 1:p.object_modes
    %                         psiq = fft2(ob{obnum}(Indy,Indx,obmode).*probes(:,:,prnum,:))/fnorm; % view in Fourier domain
    %                         chir = alpha*ifft2(p.fmask(:,:,jj).*( alpha - p.fmag(:,:,jj)./(Fq+eps) ).*psiq)*fnorm; % May not be needed for position optimization
    %                         if p.opt_flags(1) == 1,
    %                             grado{obnum}(Indy,Indx,obmode) = grado{obnum}(Indy,Indx,obmode) ...
    %                                 + sum(2*conj(probes(:,:,prnum,:)).*chir,4);
    %                         end
    %                         if p.opt_flags(2) == 1
    %                             gradp(:,:,prnum,:) = gradp(:,:,prnum,:) ...
    %                                 + 2*conj(ob{obnum}(Indy,Indx,obmode)).*chir;
    %                         end
    %                     end
    %                 end


                otherwise
                    error(['Error metric ' p.opt_errmetric 'is not defined'])
            end
        end
    end
end
func = func + initialerror;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Sieves preconditioning %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (any(smooth_gradient(:)) ~= 0)&&p.opt_flags(1)
    for obnum = 1:p.numobjs
        for obmode = 1:p.object_modes
            grado{obnum}(:,:,obmode) = conv2(grado{obnum}(:,:,obmode),smooth_gradient,'same');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Object regularization %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalized regularization to avoid the reduction of object amplitudes
% with the setting of intensity invariant
if (creg > 0)&&p.opt_flags(1)
    for obnum = 1:p.numobjs
        for obmode = 1:p.object_modes
            % Not normalized regularization
            %         func = func ...
            %             + sum(sum(   abs( ob{obnum}(2:end,1:end-1) - ob{obnum}(1:end-1,1:end-1) ).^2 ...
            %             +            abs( ob{obnum}(1:end-1,2:end) - ob{obnum}(1:end-1,1:end-1) ).^2   ));
            R = sum(sum(   abs( ob{obnum}(2:end,1:end-1,obmode) - ob{obnum}(1:end-1,1:end-1,obmode) ).^2 ...
                +          abs( ob{obnum}(1:end-1,2:end,obmode) - ob{obnum}(1:end-1,1:end-1,obmode) ).^2   ));
            norm_r = sum(sum(abs(ob{obnum}(:,:,obmode)).^2));
            func = func + creg*R/norm_r;
            if nargout > 1
                %             % Not normalized regularization
                %             grado{obnum}(2:end-1,2:end-1) = grado{obnum}(2:end-1,2:end-1) + 8*ob{obnum}(2:end-1,2:end-1) ...
                %                 -  2*ob{obnum}(1:end-2,2:end-1) - 2*ob{obnum}(3:end,2:end-1) ...
                %                 -  2*ob{obnum}(2:end-1,1:end-2) - 2*ob{obnum}(2:end-1,3:end);
                grado{obnum}(2:end-1,2:end-1,obmode) = grado{obnum}(2:end-1,2:end-1,obmode) + creg*(   (8+2*R/norm_r)*ob{obnum}(2:end-1,2:end-1,obmode) ...
                    -  2*ob{obnum}(1:end-2,2:end-1,obmode) - 2*ob{obnum}(3:end,2:end-1,obmode) ...
                    -  2*ob{obnum}(2:end-1,1:end-2,obmode) - 2*ob{obnum}(2:end-1,3:end,obmode));
            end
        end
    end
end

% normalized error, err_chi close to 1 is good result for poisson noise 
err_chi = 2*sqrt(func/prod(p.asize)/p.numpos/p.renorm^2); 

func = double(func);

if nargout > 1
    core.errorplot(err_chi);
    iteration = length(core.errorplot([]));
    verbose(2, 'Iteration # %d of %d', iteration, p.opt_iter);
    verbose(3,['Starting linesearch, Error = '  num2str(err_chi)]), 
else
    verbose(3,['Error = ' num2str(err_chi)]), 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Probe support constratint %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if p.use_probe_support&&p.opt_flags(2)
    gradp = bsxfun(@times, gradp, p.probe_mask);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Scaling preconditioning %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
avobint = 0;
if p.scale_gradient&&p.opt_flags(2)
    for ii = 1:p.numscans
        if p.share_probe
            avobint = avobint + sum( abs(grado{ii}(:)).^2 );
            if ii == p.numscans
                avobint = avobint/p.numscans;
                gradp = sqrt( avobint/sum( abs(gradp(:)).^2 ) )*gradp;
            end
        else
            gradp(:,:,ii,:) = sqrt( sum( abs(grado{ii}(:)).^2 )/sum(sum(sum( abs(gradp(:,:,ii,:)).^2 ))) )*gradp(:,:,ii,:);
        end
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Arranging gradients vector %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout > 1
    grad = []; % Optimization vector
    if p.opt_flags(1) == 1,
        for obnum = 1:p.numobjs
            grad = [grad; real(grado{obnum}(:)); imag(grado{obnum}(:))];
        end
    end
    if p.opt_flags(2) == 1,
        grad = [grad; real(gradp(:)); imag(gradp(:))];
    end
    % if flags(3) == 1,
    %     xopt = [xopt;x;y];
    % else
    %     fixed.x = x;
    %     fixed.y = y;
    % end
    if isempty(grad),
        error('At least one element of flags must be 1'),
    end
    
    %%%%%%%%%%%%%%%
    %%% Display %%%
    %%%%%%%%%%%%%%%
    
    p.error_metric.value = core.errorplot([]);
    p.error_metric.iteration = (1:size(core.errorplot([]),1));
    p.error_metric.err_metric = '-LogLik';     
    p.error_metric.method = 'ML';
    
    p.object = ob;
    p.probes = probes;
    
    if p.use_display
        if (round(mod(iteration,p.plot.interval))==0)||(iteration==1)
            p.plot.extratitlestring = sprintf(' (%dx%d) - iter %d', p.asize(2), p.asize(1), iteration);
            
            p.flat_object_used = 0;
            core.analysis.plot_results(p, 'use_display', p.use_display);
        end
    end
end

return
end

