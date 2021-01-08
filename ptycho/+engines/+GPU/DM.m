%DM Difference-Map algorithm
%
% [self, cache, psi_dash, fourier_error ] =  DM(self,par,cache,psi_dash,fourier_error, iter)
% 
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** par       structure containing parameters for the engines 
% ** psi_dash  complex projection from previous iteration, or emppty in first iteration 
% ** cache     structure with precalculated values to avoid unnecessary overhead
% ** fourier_error  array [Npos,1] containing evolution of reconstruction error 
% ** iter      current iteration number 
%
% returns:
% ++ self        self-like structure with final reconstruction
% ++ cache     structure with precalculated values to avoid unnecessary overhead
% ++ psi_dash  complex projection from previous iteration, or emppty in first iteration 
% ++ fourier_error  array [Npos,1] containing evolution of reconstruction error 
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


function [self, cache, psi_dash, fourier_error ] =  DM(self,par,cache,psi_dash,fourier_error, iter)
    import utils.*
    import math.*
    import engines.GPU.shared.*
    import engines.GPU.GPU_wrapper.*

    psi = cell(par.Nmodes, 1);
    Psi = cell(par.Nmodes, 1);

    %%%%%%%%%%%%%%%%%%%%%%%%% Difference maps algorithm %%%%%%%%%%%%%%%%%%%
    beta =  1; 
    gamma = 1; 
    relax_mask = 1; % smoothly change relaxation of the mask 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    probe_norm = norm2(self.probe{1}); 
    object_modes = length(self.object); 
    probe_amp_corr = [0,0]; 

    for ll = 1:par.object_modes
        obj_proj{ll} = Gzeros([self.Np_p, 0], true);
    end
    for ll = 1:object_modes
        obj_illum{ll,1} = Gzeros(self.Np_o);
        obj_update{ll,1} = Gzeros(self.Np_o, true);
    end
    for ll = 1:par.probe_modes
        probe_illum{ll} = Gzeros(size(self.probe{ll}));
        probe_update{ll} = Gzeros(size(self.probe{ll}), true);
    end

    % use precalculated blocks to be solved in parallel 
    indices = cache.preloaded_indices_simple{1}.indices; 
    scan_ids = cache.preloaded_indices_simple{1}.scan_ids; 
    Nind = length(indices); 
    
    p_ind = cell(max(par.probe_modes, par.object_modes),1);
    for jj = 1:Nind
        g_ind = indices{jj};         

        for ll = 1:max(par.probe_modes, par.object_modes) 
             if par.share_probe || ll > 1  % share incoherent modes 
                p_ind{ll} = 1;
             else
                 if all(scan_ids{jj} == scan_ids{jj}(1))
                    p_ind{ll} =  scan_ids{jj}(1);
                 else
                    p_ind{ll} =  scan_ids{jj};
                 end
             end
        end

        % load on GPU if needed 
        if ~par.keep_on_gpu
            for ll = 1:max(par.probe_modes, par.object_modes)
                try; psi_dash{ll,jj} = Garray(psi_dash{ll,jj}); end
            end
        end

        %% fourier propagation  
        for ll = 1:par.object_modes
            obj_proj{ll} = get_views(self.object, obj_proj{ll},1,ll, g_ind, cache, scan_ids{jj},[]);
        end
        for ll = 1:max(par.probe_modes, par.object_modes)
            probe{ll} = self.probe{min(ll,end)}(:,:,p_ind{ll});
            psi{ll} = bsxfun(@times, obj_proj{min(ll, end)},  probe{ll}); 
            if isempty(psi_dash{ll,jj}); psi_dash{ll,jj} = psi{ll}; end % initial guess 
            % P_M (P_0(psi))
            Psi{ll} = Gfun(@DM_update_psi, gamma,psi{ll}, psi_dash{ll,jj} ); 
            Psi{ll} = fwd_fourier_proj(Psi{ll}, self.modes{min(ll,end)} );  
        end  

        %% load data to GPU (if not loaded yet)
        modF = get_modulus(self, cache, g_ind,true,jj);
        mask = get_mask(self, cache.mask_indices, g_ind);

        % get intensity (modulus) on detector including different corrections
        aPsi = get_reciprocal_model(self, Psi, modF, mask,iter, g_ind, par,cache);

        if iter > 0 && (par.get_error && (mod(iter,min(20, 2^(floor(2+iter/50)))) == 0 || iter < 20) || iter == par.number_iterations )
            [fourier_error(iter,g_ind)] = get_fourier_error(modF, aPsi, [],mask);
        end

        if ~isempty(mask)
            mask = single(1)-mask;
            mask = relax_mask + (min(relax_mask, par.pfft_relaxation) - relax_mask ) * mask;
        else
            mask = par.pfft_relaxation;
        end 
        
        if  iter == 0
            % in the first iteration only find optimal scale for the probe
            probe_amp_corr(1) = probe_amp_corr(1) + Ggather(sum(modF(:).^2));
            probe_amp_corr(2) = probe_amp_corr(2) + Ggather(sum(aPsi(:).^2));
            for ll = 1:size(psi_dash,2); psi_dash{ll,jj} = []; end
            continue
        end

        %% fourier (modulus) contraint
        Psi = modulus_constraint(modF,aPsi,Psi, mask, [], par, 0 );

            
        aPsi = []; modF = []; mask = []; 


        for ll = 1:max(par.probe_modes, par.object_modes)
            Psi{ll} = back_fourier_proj(Psi{ll}, self.modes{min(end,ll)});
            psi_dash{ll,jj}  = Gfun(@DM_update, psi_dash{ll,jj} ,beta, Psi{ll},psi{ll});  
        end
        Psi = []; 

        % get back from GPU if required 
        if ~par.keep_on_gpu
            for ll = 1:max(par.probe_modes, par.object_modes)
                psi_dash{ll,jj} = Ggather(psi_dash{ll,jj});
            end
        end
    end
        
    if iter == 0
       % apply initial correction for the probe intensity and return
       probe_amp_corr = sqrt(probe_amp_corr(1) / probe_amp_corr(2)); %% calculate ratio between modF^2 and aPsi^2

       for ii = 1:par.probe_modes
           self.probe{ii} = self.probe{ii}*probe_amp_corr;
       end
       psi_dash = cell(size(psi_dash));
       verbose(2,'Probe amplitude corrected by %.3g',probe_amp_corr)
       return
    end
                  
    %% iterative solver of the overlap constraint, important for initial convergence
    for kk = 1:10
        for ll = 1:object_modes
            obj_illum{ll}(:) = 0;
            obj_update{ll}(:) = 0i;
        end
        for ll = 1:par.probe_modes
            probe_illum{ll}(:) = 0;
            probe_update{ll}(:) = 0i;
        end
         
        probe_0 = self.probe; 
        %% obj_update , obj_illum,probe_update,probe_illum is not reset to make it more stable 
        for jj = 1:length(indices)
            for ll = 1:max(par.probe_modes, par.object_modes)
                 if par.share_probe || ll > 1  % share incoherent modes 
                    p_ind{ll} = 1;
                 else
                     if all(scan_ids{jj} == scan_ids{jj}(1))
                        p_ind{ll} =  scan_ids{jj}(1);
                     else
                        p_ind{ll} =  scan_ids{jj};
                     end
                 end
            end
            for ll = 1:par.probe_modes
                probe{ll} = self.probe{ll}(:,:,p_ind{ll});
                cprobe{ll} = conj(probe{ll});
                aprobe{ll} = real(probe{ll}.*cprobe{ll});
            end
        
            % move to GPU (if not there yet)
            for ll = 1:max(par.probe_modes, par.object_modes)
                psi_dash{ll,jj} = Garray(psi_dash{ll,jj});
            end
            g_ind = indices{jj};

            for ll = 1:par.object_modes
                obj_proj{ll} = get_views(self.object, obj_proj{ll},1,ll, g_ind, cache, scan_ids{jj},[]);
            end
            if iter >= par.probe_change_start
                for ll = 1:par.probe_modes   
                    %% update probe 
                    [probe_update{ll},probe_illum{ll}] = QQ_probe(psi_dash{ll,jj}, obj_proj{min(end,ll)}, probe_update{ll},probe_illum{ll}, p_ind{ll});
                end
            end

            if iter >= par.object_change_start
                %% update object 
                [obj_update,obj_illum] = QQ_object(psi_dash{ll,jj}, obj_update,obj_illum, aprobe{min(ll,end)}, cprobe{min(ll,end)}, g_ind,min(ll, object_modes),scan_ids{jj},cache);
            end
            
            % get back from GPU if required 
            if ~par.keep_on_gpu
                for ll = 1:max(par.probe_modes, par.object_modes)
                    psi_dash{ll,jj} = Ggather(psi_dash{ll,jj});
                end
            end

        end
        for ll = 1:max([par.probe_modes, object_modes]) 
              if iter >= par.probe_change_start && ll <= par.probe_modes
                    % add some inertia to prevent oscilations 
                    self.probe{ll}  = update_probe(self, self.probe{ll} , probe_update{ll} , probe_illum{ll} , par, ll);
              end
              if iter >= par.object_change_start && ll <= object_modes
                    self.object{ll} = Gfun(@update_object, self.object{ll}, obj_update{ll}, obj_illum{ll}, cache.MAX_ILLUM(ll)*1e-4, par.probe_inertia);
              end
        end
        
        if iter > par.probe_change_start
            % if both object and probe are recontructed, solve the
            % realspace constraint iterativelly 
            min_iter = 1 + par.keep_on_gpu;  % if par.keep_on_gpu==true, then the code is so slow that it is not worthy skip the norm calculation
            min_change = 0.01; 
            if kk >= min_iter || verbose >= 2 
                dprobe = max(norm2(self.probe{1} - probe_0{1}) ./ probe_norm);
                verbose(2,'Update probe difference: %3.2f%%', dprobe*100)
            end
            if kk >= min_iter && dprobe < min_change  % change is below 1% and at least 2 iterations were done 
                break
            end   
        end
    end
    
    if verbose > 2
        for kk = 1:object_modes
            Nresid = sum2(utils.findresidues(self.object{kk}(cache.object_ROI{:})) > 0.1); 
            if Nresid > 0
                verbose(1,'Number of residua in object %i: %i', kk, Nresid)    
            end
        end
    end

end
function [probe_update,probe_illum] = QQ_probe(psi, obj_proj, probe_update,probe_illum, p_ind)
    import engines.GPU.GPU_wrapper.*
	%% update probe 
    [upd, illum] = Gfun(@QQ_probe_Gfun, psi,obj_proj);
    if size(probe_update,3)==1 % shared probe 
        probe_update = probe_update+sum(upd,3);
        probe_illum = probe_illum+sum(illum,3);
    else
        probe_update(:,:,p_ind) = probe_update(:,:,p_ind)+sum(upd,3);
        probe_illum(:,:,p_ind) = probe_illum(:,:,p_ind)+sum(illum,3);
    end
end
function norm = get_probe_norm_aux(P0, P1)
    norm = abs(P0-P1).^2; 
end

function [upd, illum] = QQ_probe_Gfun(psi,proj)
     upd = psi .* conj(proj) ;
     illum = abs(proj).^2;
end
function [obj_update,obj_illum] = QQ_object(psi,obj_update,obj_illum, aprobe, cprobe, g_ind,ll,scan_ids,cache)
    import engines.GPU.shared.*
    %% update object 
    psi = psi .* cprobe;
    [obj_update,obj_illum] = set_views_rc(obj_update,obj_illum,psi,aprobe,1,ll, g_ind, cache, scan_ids,[]);
end
function psi = DM_update_psi(gamma,psi, psi_dash )
    % real space update function for difference maps 
    psi = (1+gamma)*psi - gamma*psi_dash;
end
function psi_dash = DM_update(psi_dash,beta, psi_tmp,psi )
    % update funciton for difference maps 
    psi_dash = psi_dash  + beta * (  psi_tmp - psi  ) ;
end
function object = update_object(object, object_upd, object_illum, delta, inertia)
    % apply also some inertia in the object update 
    object = object*inertia + (1-inertia)*object_upd./ (object_illum+delta);
end
function probe = update_probe(self, probe, probe_update, probe_illum, par, probe_id)
    import engines.GPU.shared.*
    probe_new = probe_update ./  (probe_illum+1e-6); 
    % apply probe support on the first probe mode 
    if probe_id == 1
        probe_new = apply_probe_contraints(probe_new, self.modes{probe_id});
    end
    % add some inertia to prevent oscilations 
    probe = par.probe_inertia*probe + (1-par.probe_inertia)*probe_new;
end
