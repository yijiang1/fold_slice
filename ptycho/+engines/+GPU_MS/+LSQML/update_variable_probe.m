% UPDATE_VARIABLE_PROBE approximation of the OPRP method to get only the first orthogonal
% vector describing the probe 
%
% [probe, probe_evolution] = ...
%           update_variable_probe(self,probe, probe_evolution, m_probe_update,probe_update, obj_proj, chi, weights, p_ind, g_ind, cache, par)
%
% 
%
% ** self               structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** probe             [Nx,Nx,probe_modes,variable_modes] variable probe modes
% ** probe_evolution   [Npos,variable_modes] array containing evolution of the varaible modes for each position 
% ** m_probe_update    precalculated value of mean(dP,3)
% ** probe_update      probe update, ie conj(O)*chi
% ** obj_proj          [Nx,Ny,N] array, views of the object for each scan position 
% ** chi       [Nx,Ny,N] array, difference between original and updated exit-wave  
% ** weights   array of relaxation values for object pixel, reduce weight of regions with weak illumination in the variable probe calculation 
% ** p_ind      indices containg corresponding probe id for each processed position
% ** g_ind      indices corresponding to the current group that is solved in parallel 
% ** cache      structure with precalculated values to avoid unnecessary overhead
% ** par        structure containing parameters for the engines 
%
% returns:
% ++ probe             [Nx,Nx,probe_modes,variable_modes] updated variable modes 
% ++ probe_evolution   [Npos,variable_modes] updated array containing evolution of the varaible modes for each position 
%
%
% see also: engines.GPU_MS.LSQML 


%     
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
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% for LSQ-ML method 
% M. Odstrcil, A. Menzel, M.G. Sicairos,  Iterative least-squares solver for generalized maximum-likelihood ptychography, Optics Express, 2018
% for OPRP method 
%  M. Odstrcil, P. Baksh, S. A. Boden, R. Card, J. E. Chad, J. G. Frey, W. S. Brocklesby,  "Ptychographic coherent diffractive imaging with orthogonal probe relaxation." Optics express 24.8 (2016): 8360-8369
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
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



function [probe, probe_evolution] = ...
                update_variable_probe(self,probe, probe_evolution, m_probe_update,probe_update, obj_proj, chi, weights, p_ind, g_ind, cache, par)

    import math.*
    import plotting.*
    import engines.GPU_MS.GPU_wrapper.*

    
    uids = unique(p_ind);
    for kk = uids
        if length(uids) == 1  % single or shared probe between the scans 
            % avoid memory copy is possible 
            probe_update_tmp = probe_update;
            chi_tmp = chi; 
            obj_proj_tmp = obj_proj; 
            g_ind_tmp = g_ind;
            block_size = length(g_ind);
        else
            % otherwise 
            block_ind = p_ind == kk;
            probe_update_tmp = probe_update(:,:,block_ind); 
            obj_proj_tmp = obj_proj(:,:,block_ind); 
            chi_tmp = chi(:,:,block_ind); 
            g_ind_tmp = g_ind(block_ind);
            block_size= sum(block_ind);
        end


        if par.variable_probe
            % use some relaxation to avoid too faster changes 
            relax_U = min(0.1,block_size/self.Npos);
            relax_V = 1;
            % make probe_update_tmp orthogonal to the average update 
            probe_update_tmp = probe_update_tmp - m_probe_update(:,:,min(end,kk)); 
            
            for ii = 1:par.variable_probe_modes
                var_probe  = probe(:,:,kk,1+ii); 
                probe_evol = probe_evolution(g_ind_tmp,1+ii);  % evolution of the 1th SVD mode coeficient 

                [var_probe,probe_evol,probe_update_tmp] =  get_first_SVD_mode(probe_update_tmp, var_probe, probe_evol,relax_U,relax_V, obj_proj_tmp, chi_tmp);
                if ii < par.variable_probe_modes
                    % subtract projection of the updated var_probe from the
                    % probe_update_tmp to enforce orthogonality between the
                    % modes 
                    projection = sum2(probe_update_tmp .* conj(var_probe)) ./ sum2(abs(var_probe).^2); 
                    probe_update_tmp = probe_update_tmp -  projection .* var_probe; 
                end
                % return the updated vector to the probe array 
                probe(:,:,kk,1+ii) = var_probe;
                probe_evolution(g_ind_tmp,1+ii) = probe_evol;
            end

        end
        if par. variable_intensity
            % correction to account for variable intensity
            mean_probe = probe(:,:,kk,1); 
            % compare P*0 and chi to estimate best update of the intensity 
            [nom, denom] = Gfun(@get_coefs_intensity,chi_tmp, mean_probe, obj_proj_tmp);

            probe_evolution(g_ind_tmp,1)  = probe_evolution(g_ind_tmp,1) + 0.1* squeeze(Ggather(sum2(nom)./ sum2(denom)));
        end
    end

    
    if any(g_ind==1) && utils.verbose() > 3
        self.probe{1} = probe; 
        self.probe_evolution = probe_evolution; 
        plot_variable_probe(self, par)
        drawnow
    end  
end

function [var_probe,probe_evol, probe_update] = get_first_SVD_mode(probe_update, var_probe, probe_evol, relax_U,relax_V, obj_proj, chi)
    import engines.GPU_MS.GPU_wrapper.*
    import math.*
    import plotting.*
    


    % get a weighting function => avoid effect of too strong noise
    % around edges of the reconstructed region => improve robustness againts outliers 
%             weights = weights / max2(weights);
%             weight_proj = get_views( weights,  Gzeros(size(chi)),1,1, g_ind, cache);
    weight_proj = 1;

    %% calculate terms needed to calculate update of the variable probe
    % => U term in SVD decomposition 
    [resid, proj, probe_update] = Gfun(@get_SVD_update,probe_update, weight_proj, var_probe, reshape(probe_evol,1,1,[]), norm(probe_evol));
    % get update the variable probe 
    var_probe_upd = mean( resid .* mean2(proj), 3);
    % apply update, prevent too large changes at the beginning of the covergence 
    var_probe = var_probe + relax_U*var_probe_upd / norm2(var_probe_upd);

    %% equivalent but much slower code 
    %[U,S,V] = svd(reshape(weight_proj.*(probe_update - m_probe_update), prod(self.Np_p),[]), 0);
    %var_probe = var_probe + relax_U*reshape(U(:,1), self.Np_p); 

    % keep the eigenprobe normalized 
    var_probe = var_probe ./ norm2(var_probe);

    %% calculate optimal OPRP evolution coeficients 
    [num, denum] = Gfun(@get_SVD_evol,var_probe, obj_proj, chi);
    num =Ggather(mean2(num)); 
    denum = Ggather(mean2(denum)); 
    %  perform relaxed update => improve robustness againts outliers 
    probe_evol_upd = squeeze(num  ./ (denum + 0.1*mean(denum,3))); % caclulate regularized update 
    % add to the coefficients that are already used in the currently
    % used probe{1} variable 

    probe_evol = probe_evol + relax_V*probe_evol_upd; 


            
end


%% GPU kernel merging

% SVD approximation => calculation of U 
function [resid, proj, probe_update] = get_SVD_update(probe_update, weight_proj, var_probe, probe_evol, probe_evol_norm)
    resid = weight_proj .* probe_update; 
    proj =  (real(conj(resid) .* var_probe)+ probe_evol) / probe_evol_norm^2;
end

% SVD approximation => calculation of S*V
function [num, denum] = get_SVD_evol(var_probe, obj_proj, chi)
    psi = var_probe .* obj_proj;
    denum = abs( psi ).^2; 
    num =  real(chi .* conj(psi));  
end

function [nom1, denom1] = get_coefs_intensity(xi, P, O)
    OP = O.*P;
    nom1 = real(conj(OP) .* xi);
    denom1 = abs(OP).^2; 
end


