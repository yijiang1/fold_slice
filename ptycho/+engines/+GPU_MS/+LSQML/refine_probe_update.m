% REFINE_PROBE_UPDATE calculate improved update direction
% apply "overlap" constraint to get better estimate  of the update direction
% also updates the variable probe estimate 
%
% [ self,m_probe_update, probe_update, cache]  = ...
%                 refine_probe_update(self, obj_proj, probe_update, chi,layer_ids,probe_id,p_ind,g_ind, par, cache)
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** obj_proj  [Nx, Ny, N] array, estimate of the object update for each scan position, ie conj(P)*chi
% ** probe_update      cell of object sized arrays containg previous optimal updates 
% ** chi      indices corresponding to the current group that is solved in parallel 
% ** layer_ids  id of the solved layer for multilayer ptycho 
% ** probe_id   id of incoherent probe mode 
% ** p_ind      indices containg corresponding probe id for each processed position
% ** g_ind      indices corresponding to the current group that is solved in parallel 
% ** par       structure containing parameters for the engines 
% ** cache     structure with precalculated values to avoid unnecessary overhead
%
% returns:
% ++ object_upd_sum      cell of object sized arrays containg updated optimal update
% ++ object_update_proj  [Nx, Ny, N] array, estimate of the refiend object update for each scan position,
% ++ cache                structure with precalculated values to avoid unnecessary overhead
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



function [ self,m_probe_update, probe_update, cache]  = ...
                refine_probe_update(self, obj_proj, probe_update, chi,layer_ids,probe_id,p_ind,g_ind, par, cache)
    % get probe update direction 
    import engines.GPU_MS.shared.*
    import engines.GPU_MS.LSQML.*
    import engines.GPU_MS.GPU_wrapper.*
    import math.*
    import utils.*

    if layer_ids > 1 % in case of multilayer object 
        m_probe_update = []; 
        return
    end
    
    if (probe_id == 1 && par.apply_subpix_shift && isinf(self.z_distance(end))) || is_used(par,'fly_scan')
        probe_update  = apply_subpx_shift(probe_update , -self.modes{min(end,probe_id)}.sub_px_shift(g_ind,:) );
    end
    if probe_id == 1
        probe_update = apply_subpx_shift_fft(probe_update, -self.modes{min(end,probe_id)}.probe_fourier_shift(g_ind,:));
    end

  
    if par.share_probe || length(unique(p_ind)) == 1 
       % BETTER WAY: assume that sum(|obj_proj|^2,3) is close to 1 
       % and additionally use weighting based on confidence given by illum_sum_0
        % => apriory weighting giving less importance to the less
        % illuminated regions
        %weight_proj = cache.illum_sum_0{1} ./ (cache.illum_sum_0{1}+0.01*cache.MAX_ILLUM(1));
        %weight_proj = get_views({weight_proj},[],1,1, g_ind, cache);
        %m_probe_update = mean( weight_proj.* probe_update,3);
        % or originally was used simple average , good for object >> probe 
        m_probe_update = mean(probe_update,3);  % calculate single update for all current positions
    else  % unshared probe and multiple scans in one group (ie shared object)
        probe_ids = unique(p_ind); 
        m_probe_update = Gzeros([self.Np_p,length(probe_ids)], true);
        for probe_id =probe_ids(:)'
            m_probe_update(:,:,probe_id) = mean(probe_update(:,:,p_ind == probe_id),3); % calculate one update for each scan
        end
    end

    if (par.variable_probe || par.variable_intensity) && probe_id == 1  
        % ORTHOGONAL PROBE RELAXATION (OPRP) EXTENSION - allow
        % variable probe wavefront 
        %  Odstrcil, M., et al. "Ptychographic coherent diffractive imaging with orthogonal probe relaxation." Optics express 24.8 (2016): 8360-8369.
        % iterate over all sub probes 
        [self.probe{probe_id}, self.probe_evolution] = ...
            update_variable_probe(self, self.probe{probe_id}, self.probe_evolution, m_probe_update,probe_update, obj_proj, chi,cache.illum_sum_0{probe_id}, p_ind, g_ind, cache, par);
    end

%     % apply probe constraints 
%     if probe_id == 1 && (check_option(self,'probe_support') || check_option(self,'probe_support_fft'))
%         m_probe_update = apply_probe_contraints(m_probe_update, self.modes{probe_id});
%     end

end