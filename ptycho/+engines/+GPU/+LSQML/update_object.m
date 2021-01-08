% UPDATE_OBJECT calculate improved update direction
% apply "overlap" constraint to get better estimate  of the update direction
% 
% [object, object_upd_sum] = update_object(self, object, object_upd_sum, layer,object_ids, g_ind, scan_ids, par, cache, beta_object)
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** object    cell of object arrays
% ** object_upd_sum      cell of object sized arrays containg previous optimal updates 
% ** layer_ids  id of the solved layer for multilayer ptycho 
% ** object_ids id of incoherent object mode (not implemented)
% ** g_ind      indices corresponding to the current group that is solved in parallel 
% ** scan_ids   determines to which scan correponds each of the position 
% ** par       structure containing parameters for the engines 
% ** cache     structure with precalculated values to avoid unnecessary overhead
% ** beta_object  (scalar) relaxation parameter of the update step 
%
% returns:
% ++ object    cell of object arrays, after update 
% ++ object_upd_sum      cell of object sized arrays containg updated optimal update
%
%
% see also: engines.GPU.LSQML 


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
% 5.	LICENSEE agrees that it shall make the foobject_idswing acknowledgement in any publication resulting from the use of the PROGRAM or any translation of the code into 
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
%       make any copies except for its internal use, without prior written consent of LICENSOR. LICENSEE agrees to place the foobject_idswing copyright notice on any such copies: 
%       © All rights reserved. PAUL SCHERRER INSTITUT, Switzerland, Laboratory for Macromolecules and Bioimaging, 2017. 
% 8.	The LICENSE shall not be construed to confer any rights upon LICENSEE by implication or otherwise except as specifically set forth herein.
% 9.	DISCLAIMER: LICENSEE shall be aware that Phase Focus Limited of Sheffield, UK has an international portfolio of patents and pending applications which relate 
%       to ptychography and that the PROGRAM may be capable of being used in circumstances which may fall within the claims of one or more of the Phase Focus patents, 
%       in particular of patent with international application number PCT/GB2005/001464. The LICENSOR explicitly declares not to indemnify the users of the software 
%       in case Phase Focus or any other third party will open a legal action against the LICENSEE due to the use of the program.
% 10.	This Agreement shall be governed by the material laws of Switzerland and any dispute arising out of this Agreement or use of the PROGRAM shall be brought before 
%       the courts of Zürich, Switzerland. 



function [object, object_upd_sum] = update_object(self, object, object_upd_sum, layer,object_ids, g_ind, scan_ids, par, cache, beta_object)
    import engines.GPU.shared.*
    import engines.GPU.GPU_wrapper.*
    import math.*
    import utils.*
    
    if par.share_object
        obj_ids = 1;
    else
        obj_ids = unique([scan_ids{:}]);
    end
    % take single optimal value, works well for most of samples 
    % it is possible to use different weighting for each scan position but
    % it may become less stable in some cases -> robusness is preferred 
    
    % in case of the MLc method take minimum of the LSQ updates from all
    % subsets 
    
    if is_method(par, 'MLc')
        % preconditioner should be applied on the total sum of all object_upd_sum
        if par.delta_p > 0 % && par.Nlayers == 1
            for kk = obj_ids
                object_upd_sum{kk,layer} = Gfun(@object_sum_update_Gfun, object_upd_sum{kk,layer}, cache.illum_sum_0{kk},cache.MAX_ILLUM(kk)*(par.delta_p));
            end
        end
    end
    
    % calculate optimal step, apply at least different step for each object (scan)
    for i = 1:length(g_ind)
        for kk = obj_ids
            ind = g_ind{i}(scan_ids{i} == kk); 
            if isempty(ind)
                beta_object_avg(i,kk) = nan;
            else
                beta_object_avg(i,kk) = trimmean(beta_object(ind), 10);
            end
        end
    end
       
    % take the most pesimistic estimate of the per object 
    beta_object_avg=nanmin(beta_object_avg,[],1);

    % update each of the objects separately 
    for kk = obj_ids
        if beta_object_avg(kk) > 0
            object_upd_sum{kk,layer} = object_upd_sum{kk,layer}*beta_object_avg(kk); 
            object{kk,layer} =  object{kk,layer}+object_upd_sum{kk,layer}; 
        end
    end
    
    if verbose()> 3
        % show applied subsets (update amplitude) and probe update amplitude 
        plotting.smart_figure(11)
        Nobj = size(object_upd_sum,1); 
        for ll = 1:Nobj
            subplot(Nobj,2,1+Nobj*(ll-1))
            cla()
            o = object_upd_sum{ll,layer}(cache.object_ROI{:}); 
            o = min(abs(o), quantile(abs(o(:)), 0.999)) .* o ./ abs(o); 
            plotting.imagesc3D(o);
            axis image xy off 
            hold all
            fprintf('Object update norm: %g\n', norm2(object_upd_sum{ll,layer}(cache.object_ROI{:})))
%             try
%                 for k = 1:length(g_ind)
%                     for i = unique(scan_ids{k})
%                         plot(self.probe_positions_0(g_ind{k}(scan_ids{k}==i),1)+self.Np_o(2)/2,self.probe_positions_0(g_ind{k}(scan_ids{k}==i),2)+self.Np_o(1)/2, '.')
%                     end
%                 end
%             end
        end
        title('Object update')
    end
            
    
end

function object_upd_sum = object_sum_update_Gfun(object_upd_sum, obj_illum_sq_sum, max)
    % final update is just weighted mean of the updared object views  
    object_upd_sum =  object_upd_sum ./ sqrt(obj_illum_sq_sum.^2+ max.^2);
end
