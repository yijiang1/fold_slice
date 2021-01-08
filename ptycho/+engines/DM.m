%DM Difference-Map algorithm
%
%   Publications most relevant to the Difference-Map implementation
%       + P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, 
%       "High-Resolution Scanning X-ray Diffraction Microscopy," Science 321, 379-382 (2008)
%       + P. Thibault, M. Dierolf, O. Bunk, A. Menzel, F. Pfeiffer,
%       "Probe retrieval in ptychographic coherent diffractive imaging,"
%       Ultramicroscopy 109, 338–343 (2009)

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

function [ p, fdb ] = DM( p )

import math.*
import utils.verbose
import utils.pshift


global proj1_time  objproj_time probeproj_time elsewheretime proj2_time plot_time

fdb.status = core.engine_status;


% mex or matlab - object_update / probe_update / Fourier_loop
if ~isfield(p, 'use_mex')
    p.use_mex = zeros(1,3);
elseif numel(p.use_mex) == 1
    if any(p.use_mex)
        verbose(3, 'Using mex files for DM.')
    end
    p.use_mex = repmat(p.use_mex,1,3);
end


% Starting iterate

% iter is a large array and is only needed for Matlab Difference map
iter = cell(p.numscans,1);
fmag = cell(p.numscans,1);
for ii = 1:p.numscans
    iter{ii} = complex(zeros([p.asize, length(p.scanidxs{ii}), p.probe_modes, p.object_modes]));
    fmag{ii} = double(p.fmag(:,:,p.scanidxs{ii})); % can cause memory duplication if it was not double, migrate to single in future !!!
end

if p.object_modes == 1 && p.numprobs == 1
    obj_proj = complex(zeros([p.asize, length(p.scanidxs{1})]));
end

if  p.probe_mask_bool
    
    if p.probe_mask_use_auto
        verbose(3, 'Using a probe mask from probe autocorrelation.');
        to_threshold = -real(auto);
    else
        verbose(3, 'Using a circular probe mask.');
        [x,y] = meshgrid(-p.asize(2)/2:floor((p.asize(2)-1)/2),-p.asize(1)/2:floor((p.asize(1)-1)/2));
        to_threshold = (x.^2 + y.^2);
        clear x y
    end

    p.probe_mask  = to_threshold < quantile(to_threshold(:), p.probe_mask_area );
else
    p.probe_mask = 1;
end
    

% check what is the actual size of the object 
for i = 1:p.numobjs
    p.object_size(i,:) = size(p.object{i});
    ob{i} = double(p.object{i});
end
% move to doubles (not needed)
p.probes = double(p.probes);
p.fmag = double(p.fmag); 


for obnum = 1:p.numobjs
    avob{obnum} = zeros([p.object_size(obnum,:) p.object_modes]);
end

% get views from probes and objects
for ii = 1:p.numscans
    prnum = p.share_probe_ID(ii); 
    obnum = p.share_object_ID(ii); 
                            
    if p.object_modes == 1 && p.numprobs == 1
        % faster version without extra memory allocation 
        obj_proj = core.get_projections(p, ob{obnum}, ii,obj_proj);
        iter{ii} = bsxfun(@times, p.probes, obj_proj);
    else
        % general version 
        for obmode = 1:p.object_modes
            obj_proj = core.get_projections(p, ob{obnum}(:,:,obmode), ii);
            iter{ii}(:,:,:,:,obmode) = bsxfun(@times, p.probes(:,:,prnum,:), obj_proj);
        end
    end
end

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
        p.userflatind = (p.object_flat_region == 1);% Find indices of flat region
    end
else
    p.userflatregion = false;
end

for obnum = 1:p.numobjs
    avob{obnum} = zeros([p.object_size(obnum,:) p.object_modes]);
end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Main Difference map loop %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    
    % A power bound (relaxed Fourier) that scales with number of photons per diffraction pattern
    p.power_bound = p.count_bound*p.renorm^2;
   
    
    
        % Prepare statistics
    proj1_time = 0;
    proj2_time = 0;
    plot_time = 0;
    objproj_time =   0;
    probeproj_time = 0;
    elsewheretime =  0;

    
    cfact_temp = p.probe_regularization *p.numpts;
    if p.share_probe
        cfact = zeros(1,p.numprobs);
        for ii=1:p.numscans
            cfact(p.share_probe_ID(ii)) = cfact(p.share_probe_ID(ii)) + cfact_temp(ii);
        end
    else
        cfact = cfact_temp;
    end
    
    err = 0;
    rfact = 0;
    numav = 0;
    
    for it=1:p.number_iterations
        verbose(2, 'Iteration # %d of %d',it, p.number_iterations);
        

        
        
        % 1. Overlap projection - a large loop where probe and object are refined.
        verbose(3, ' - projection 1: overlap constraint - ');
        proj1tic = tic;
        
        % Check once each iteration to see whether the named file exists,
        % signalling a user break.
        fid = fopen(p.io.break_check_name);
        if (fid ~= -1)
            fprintf('The file %s exists.',p.io.break_check_name);
            fprintf('It will be deleted now and then the program terminates.\n');
            delete(p.io.break_check_name);
            sprintf('Terminating program upon user request.\n');
            break
        end
        
        % The simple iterative scheme
        prch0 = 0;
        breakprobeloop = 0;
        for inner=1:10
            if ~breakprobeloop
                cprobes = conj(p.probes);
                for obnum = 1:p.numobjs
                    pr_nrm{obnum} = 1e-8 * ones(p.object_size(obnum,:));
                    %                         % This could be zero, but better make it a very small number instead to avoid eventual divisions by 0.
                    ob{obnum} = 1e-8*(1+1i)*ones([p.object_size(obnum,:)  p.object_modes]);
                    %pr_nrm{obnum} = 1e-0*ob{obnum};
                    % This could be zero, but better make it a very small number instead to avoid eventual divisions by 0.
                    %ob{obnum} = 1e-0*ob{obnum};
                end
                for ii = 1:p.numscans
                    objtic = tic;
                    prnum = p.share_probe_ID(ii); 
                    obnum = p.share_object_ID(ii); 
                    % Decide which matlab code to use
                    if (p.probe_modes == 1)&&(p.object_modes == 1)&&p.use_mex(1)
                        % Mex code object loop
                        engines.DM.object_update_norm(iter{ii},...
                            cprobes(:,:,prnum),ob{obnum},pr_nrm{obnum},int32(p.positions(p.scanidxs{ii},:)),...
                            int32(p.numpts(ii)));
                    else
                        % Pure matlab version
                        if p.object_modes == 1 && p.numprobs == 1
                            % faster version without memory allocation 
                            cprobe = cprobes; 
                            obj_update = sum(bsxfun(@times, cprobe, iter{ii}),4);
                            ob{obnum} = core.set_projections(p, ob{obnum}, obj_update , ii);
                        else
                            cprobe = cprobes(:,:,prnum,:); 
                            for obmode = 1:p.object_modes
                                obj_update = sum(bsxfun(@times, cprobe, iter{ii}(:,:,:,:,obmode)),4);
                                ob{obnum}(:,:,obmode) = core.set_projections(p, ob{obnum}(:,:,obmode), obj_update , ii);
                            end 
                        end
                        pr_nrm{obnum} = core.set_projections(p, pr_nrm{obnum}, sum(abs(cprobe).^2,4) , ii);
                    end
                    objproj_time = objproj_time + toc(objtic);
                end               
                                
                for obnum = 1:p.numobjs
                    ob{obnum} = bsxfun(@rdivide, ob{obnum}, pr_nrm{obnum});
                                       
                    if p.userflatregion % Not supported with modes
                        ob{obnum}(userflatind) = mean(ob{obnum}(userflatind));
                    end
                    
                    elsewheretic = tic;
                    if p.clip_object
                        aob = abs(ob{obnum});
                        too_high = (aob > p.clip_max);
                        too_low = (aob < p.clip_min);
                        ob{obnum} = (1-too_high).*(1-too_low).*ob{obnum} + (too_high.*p.clip_max + too_low*p.clip_min).*ob{obnum}./aob;
                    end
                    elsewheretime = elsewheretime + toc(elsewheretic);
                    
                end
                if p.probe_change_start >= it
                    breakprobeloop = 1;
                else
                    % Defining the new probes (regularization)
                    nprobes = bsxfun(@times,p.probes, reshape(cfact,1,1,[]));
                    pr_denoms = bsxfun(@times,ones(p.asize), reshape(cfact,1,1,[]));
                    for ii = 1:p.numscans
                        prnum = p.share_probe_ID(ii); 
                        obnum = p.share_object_ID(ii); 
                        
                        probetic = tic;
                        nprobe =     nprobes(:,:,prnum,:); % Current scan probe, third index is modes
                        pr_denom = pr_denoms(:,:,prnum); % Current scan denominator, no mode index
                        % Decide which code to use
                        if (p.probe_modes == 1)&&(p.object_modes == 1)&&p.use_mex(2)
                           % Mex code
                            engines.DM.probe_update_norm(iter{ii},...
                               squeeze(nprobe),ob{obnum},pr_denom,...
                                int32(p.positions(p.scanidxs{ii},:)),int32(p.numpts(ii)));
                        else
                            % Pure matlab code
                            if p.object_modes == 1 &&  p.numprobs == 1
                                % faster version without memory allocation
                                obj_proj = core.get_projections(p, ob{obnum}, ii, obj_proj);
                                nprobe = nprobe + sum(bsxfun(@times,iter{ii}, conj(obj_proj)), 3);% sum over positions 
                                pr_denom = pr_denom + sum(abs(obj_proj).^2,3);% sum over positions 
                           else
                               for obmode = 1:p.object_modes
                                    obj_proj = core.get_projections(p, ob{obnum}(:,:,obmode), ii);
                                    nprobe = nprobe + sum(bsxfun(@times,iter{ii}(:,:,:,:,obmode), conj(obj_proj)),3);  % sum over positions 
                                    pr_denom = pr_denom + sum(abs(obj_proj).^2,3); % sum over positions 
                                end
                           end
                        end
                        probeproj_time = probeproj_time + toc(probetic);

                        nprobes(:,:,prnum,:) = nprobe;
                        pr_denoms(:,:,prnum) = pr_denom; 
                    end
                    
                    probe_new = bsxfun(@rdivide, nprobes, pr_denoms);
                    probe_new = bsxfun(@times, p.probe_mask , probe_new); 
                                        
                    % get relative residuum 
                    prch = squeeze(sqrt(sum(sum(sum(abs(p.probes - probe_new).^2,1),2),4) ./ sum(sum(sum(abs(p.probes).^2,1),2),4))); 
                    for prnum = 1:p.numprobs
                        verbose(3, 'Change in probe %d: %3.2g%%',prnum,prch(prnum)*100);
                    end
                    p.probes = probe_new; 

                    if all(prch < 0.01)
                        breakprobeloop = 1;
                    end
                end
            end
        end
        
        proj1_time = proj1_time + toc(proj1tic);
        
        er2 = 0;
        % 2. Fourier projection + complete difmap loop
        verbose(3, ' - projection 2: Fourier modulus constraint - ');
        

        
        
        % perform normalization of probe (avoid probe - object scaling ambiguity)
        if (it < p.average_start && p.remove_scaling_ambiguity)
            pnorm = math.norm2(p.probes);

            if length(unique(p.share_probe_ID)) == p.numscans && ...
               length(unique(p.share_object_ID)) == p.numscans

                p.probes = p.probes ./ pnorm; 
                for ii = 1:p.numscans
                    ob{ii} = ob{ii} .* pnorm(ii); 
                end
            else
                pnorm = mean(pnorm); 
                p.probes = p.probes / pnorm; 
                for ii = 1:p.numscans
                    ob{ii} = ob{ii} * pnorm; 
                end
            end
        end
        

        tic;
        if (p.probe_modes == 1)&&(p.object_modes == 1)&&p.use_mex(3)
            for ii = 1:p.numscans
                prnum = p.share_probe_ID(ii); 
                obnum = p.share_object_ID(ii); 
                
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
        
                fmask_per_scan = ndims(p.fmask) == 3; 
                % Mex code
                %         Fourier_DM_loop_par2(iter,probe,ob{obnum},p1,p2,f,ph,df,double(fmask),fmag,double(p.power_bound),er2,rf,rf_nrm,af, fdev, fdev2, fmaski,int32(positions),int32(sum(p.numpts)),int32(fmask_per_scan),int32(compute_rfact));               
                engines.DM.Fourier_DM_loop_par2(iter{ii},...
                    (p.probes(:,:,prnum)),(ob{obnum}),p1,p2,f,ph,df,...
                    double(p.fmask(:,:,p.scanidxs{ii})),...
                    fmag{ii},double(p.power_bound),...
                    er2,rf,rf_nrm,af, fdev, fdev2, fmaski,int32(p.positions(p.scanidxs{ii},:)),...
                    int32(p.numpts(ii)),int32(fmask_per_scan),int32(p.compute_rfact));


                % What the Mex function does:
                %                 fnorm = sqrt(a2);
                %                 rf = 0;
                %                 rf_nrm = 0;
                %                 if ~p.fmask_per_scan
                %                     fmaski = fmask;
                %                 end
                %                 for jj=p.scanidxs{ii}
                %                     if fmask_per_scan
                %                         fmaski = fmask(:,:,jj);
                %                     end
                %                     Indy = p.positions(jj,1) + (1:p.asize(1));
                %                     Indx = p.positions(jj,2) + (1:p.asize(2));
                %                     p1 = p.probes(:,:,prnum) .* ob{obnum}(Indy, Indx);
                %
                %                     f = fft2( 2*p1 - iter(:,:,jj) )/fnorm;
                %                     af = abs(f);
                %                     ph = f ./ (af+1e-10);
                %                     fdev = af - fmag(:,:,jj);
                %                     fdev2 = fmaski.*fdev.^2;
                %                     power = sum(sum(fdev2))/a2;
                %                     if power > p.power_bound
                %                         renorm = sqrt(p.power_bound / power);
                %                         af = af.*(1-fmaski) + fmaski.*(fmag(:,:,jj) + fdev * renorm);
                %                     end
                %                     p2 = fnorm*ifft2(af .* ph);
                %
                %                     df = p2 - p1;
                %                     iter(:,:,jj) = iter(:,:,jj) + df;
                %
                %                     er2 = er2 + sum(sum(abs(df).^2));
                %
                %                     if p.compute_rfact
                %                         rf = rf + sum(sum(abs(abs(fft2(p1))/fnorm - fmag(:,:,jj))));
                %                         rf_nrm = rf_nrm + sum(sum(fmag(:,:,jj)));
                %                     end
                %                 end
            end
        else
            % Pure Matlab version
%             verbose(2,'Matlab Fourier loop')
            % The solution seems equivalent down to the displayable
            % digits, but when starting optimization the error metric
            % of the current guess is slightly different
 
            for ii = 1:p.numscans
                prnum = p.share_probe_ID(ii); 
                obnum = p.share_object_ID(ii); 
                    
                fnorm = sqrt(prod(p.asize));

                rf = 0;
                rf_nrm = 0;
                                 

                p1 = zeros([p.asize,length(p.scanidxs{ii}),p.probe_modes,p.object_modes]);
                if p.object_modes == 1 && p.numprobs == 1
                    % faster version without memory allocation
                    obj_proj = core.get_projections(p, ob{obnum}, ii,obj_proj);
                    p1 = bsxfun(@times,p.probes,obj_proj);
                else
                    for obmode = 1:p.object_modes
                        obj_proj = core.get_projections(p, ob{obnum}(:,:,obmode), ii);
                        p1(:,:,:,:,obmode) = bsxfun(@times,p.probes(:,:,prnum,:),obj_proj);
                    end
                end
                f = fft2( 2*p1 - iter{ii} )/fnorm;
                af = abs(f); % Amplitude of f before projection
                ph = f ./ (af+1e-3);
                fmag_target = bsxfun(@times, af, fmag{ii}./sqrt(sum(af.^2,4)));
                fdev = af - fmag_target;
                if size(p.fmask,3)==p.numpos
                    fmaski = p.fmask(:,:,p.scanidxs{ii}); 
                else
                    fmaski = p.fmask;
                end
                af = bsxfun(@times, af, 1-fmaski) + bsxfun(@times, fmaski, fmag_target + fdev * p.pfft_relaxation);
                p2 = fnorm*ifft2(af .* ph);
                clear ph af 
                df = p2 - p1;
                clear p1 p2
                iter{ii} = iter{ii} + df;
                er2 = sum2(squeeze(sum2(abs(df).^2)));  % faster than (:)
                clear df

            end
        end
        
        if p.center_probe
            % Find probe centroid
            for ii = 1:p.numscans
                prnum = p.share_probe_ID(ii); 
                obnum = p.share_object_ID(ii);

                
                Iprobe = abs(p.probes(:,:,prnum,:)).^2;
                [probe_c2, probe_c1]=center(Iprobe);

                % shift all arrays accordingly
                if (probe_c1 ~= 1) || (probe_c2 ~= 1)
                    verbose(3,'Shifting all arrays by (%d, %d)', probe_c1,probe_c2);
                    p.probes(:,:,prnum) = pshift(p.probes(:,:,prnum,:),[probe_c1 probe_c2]);
                    ob{obnum} = pshift(ob{obnum},[probe_c1 probe_c2]);
                    avob{obnum} = pshift(avob{obnum},[probe_c1 probe_c2]);
                    for i=p.scanidxs{ii}
                        iter(:,:,i,:) = pshift(iter(:,:,i,:), [probe_c1 probe_c2]);
                    end
                end
            end
        end
        
        proj2_time = proj2_time + toc();
        err(it) = 2*sqrt(er2/(prod(p.asize)*sum(p.numpts)));
        if p.compute_rfact
            rfact(it) = rf/rf_nrm;
            verbose(3, 'Error: %12.3f \t R-factor: %12.3f %%',err(it),100*rfact(it));
        else
            verbose(3,'Error: %12.3f',err(it));
        end
        
        if (it >= p.average_start) && mod(it, p.average_interval)==0
            for obnum = 1:p.numobjs
                avob{obnum} = avob{obnum} + ob{obnum};
            end
            numav = numav + 1;
        end
        
        
        p.error_metric.iteration = 1:it ;
        p.error_metric.value = err(1:it);
        p.error_metric.err_metric = 'RMS';
        p.error_metric.method = p.name;
        p.object = ob;

        tic;
        if (round(mod(it,p.plot.interval))==0)||(it==1) && p.use_display
            p.plot.extratitlestring = sprintf(' (%dx%d) - iter %d', p.asize(2), p.asize(1), it);
            core.analysis.plot_results(p);
        end
        plot_time = plot_time + toc();
        
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% end main difference map loop %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    verbose(3, 'Finished difference map');
    verbose(3, 'Time elapsed in projection 1: %f seconds', proj1_time);
    verbose(3, '             in object projection: %f seconds', objproj_time);
    verbose(3, '             in probe projection: %f seconds', probeproj_time);
    verbose(3, 'Time elapsed in projection 2: %f seconds', proj2_time);
    verbose(3, 'Time spent plotting: %f seconds', plot_time);
        
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
        prnum = p.share_probe_ID(ii); 
        obnum = p.share_object_ID(ii);
        for i=p.scanidxs{ii}
            Indy = round(p.positions(i,1)) + (1:p.asize(1));
            Indx = round(p.positions(i,2)) + (1:p.asize(2));
            fcalc = abs(fftn(avob{obnum}(Indy,Indx).*p.probes(:,:,prnum)))/sqrt(prod(p.asize));
            av_rfact = av_rfact + sum(sum(abs(fcalc - p.fmag(:,:,i))));
            av_rfact_nrm = av_rfact_nrm + sum(sum(p.fmag(:,:,i)));
        end
    end
    av_rfact = av_rfact / av_rfact_nrm;
    
    p.object = avob;

    
    % save additional feedback into engine's fdb variable
    fdb.rfact = rfact;
    fdb.av_rfact = av_rfact;




end
