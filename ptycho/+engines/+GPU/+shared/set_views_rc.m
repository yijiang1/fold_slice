% SET_VIEWS_RC Reduce stack of projections into one shared object
% ! update complex and real object at once 
%
% [obj_update,obj_illum] = set_views_rc(obj_update,obj_illum, psi,aprobe2,layer,object_id, indices, cache, scan_ids, skip_ind, object_modes)
%
% ** obj_update   [Nx_o, Ny_o] array or cells containing object 
% ** obj_illum    [Nx_o, Ny_o] array or cells containing illumination sum  
% ** obj_proj   [Nx_p, Ny_p, N] preallocated array for the views
% ** psi        complex valued exit wave (psi ~ P*O)
% ** aprobe2     illumination intensity patch 
% ** layer_ids  id of the solved layer for multilayer ptycho 
% ** object_id  id of the object, ie scan or incoherent mode 
% ** indices   processed positions 
% ** cache     structure with precalculated values to avoid unnecessary overhead
% ** scan_ids   determines to which scan correponds each of the position 
% ** skip_ind   list of indices to be skipped 
% ** object_modes   number of incoherent object modes 
%
% returns:
% ++ object    reduced sum of the views




%     

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



function [obj_update,obj_illum] = set_views_rc(obj_update,obj_illum, psi,aprobe2,layer,object_id, indices, cache, scan_ids, skip_ind)

    import engines.GPU.GPU_wrapper.*
    import engines.GPU.shared.*
    import utils.verbose
    global use_gpu 
        
   if nargin < 10
        skip_ind = [];
    end
    if nargin > 8 && ~isempty(scan_ids)
       
        %% !!! call recursivelly -> wrapper for multiscan version !!!!
        
        % get unique IDs of the scans
        if isempty(scan_ids)
            unq_scans = []; 
        elseif all(scan_ids == scan_ids(1))
            unq_scans = scan_ids(1);
        else
            unq_scans = unique(scan_ids);
        end
        
        if length(unq_scans)> 1 || (iscell(obj_update) &&  length(obj_update) > 1)
          if ~isempty(use_gpu) && use_gpu && isa(obj_update{1}, 'gpuArray') && isa(psi, 'gpuArray') 
               if size(obj_update,1) == 1
                    % shared object or single object 
                    ind_ok{1} = uint16(1:length(indices));
               else
                    for kk = 1:size(obj_update,1)
                        ind_ok{kk} = uint16(find(scan_ids == kk));
                    end
               end
               % feed data directly to the GPU mex without splitting 
               [obj_update(:,layer),obj_illum(:,layer)] = set_views_gpu_rc(obj_update(:,layer),obj_illum(:,layer),psi, aprobe2,cache.oROI_s{min(object_id,end)},indices, ind_ok );
           else
                % if GPU not available use this "wrapper" around single
                % set_projection function 
                for kk = unq_scans
                    ind = scan_ids == kk;
                    skip_ind = indices(~ind);  % avoid going through these indices 
                    [obj_update{kk,layer},obj_illum{kk,layer}] = set_views_rc(obj_update{kk,layer},obj_illum{kk,layer},psi,aprobe2, layer, object_id, indices, cache,[],skip_ind);
                end
          end
          return 
        end
    end
     
    
    is_cell =  iscell(obj_update);
    if is_cell
        obj_update_0 = obj_update; 
        obj_illum_0 = obj_illum;
        obj_update = obj_update{object_id, min(end,layer)};
        obj_illum = obj_illum{object_id, min(end,layer)};
    end
    
    if ~isempty(cache.skip_ind) && ~isempty(skip_ind)
        ind_ok = uint16(find(~ismember(indices, [cache.skip_ind,skip_ind])));  % skip wrong patterns 
    else
        ind_ok = uint16(1:length(indices)); 
    end
    
    if ~isempty(use_gpu) && use_gpu && isa(obj_update, 'gpuArray') && isa(psi, 'gpuArray')
        %% USE CUDA MEX FOR GPU 
        [obj_update,obj_illum] = set_views_gpu_rc(obj_update,obj_illum, psi,aprobe2,cache.oROI_s{min(object_id,end)},indices, ind_ok); 
    else  
        %% USE CPU 
        obj_illum = set_views(obj_illum , aprobe2, layer, object_id,indices, cache,scan_ids, skip_ind);
        obj_update = set_views(complex(obj_update),complex(psi), layer, object_id, indices, cache,scan_ids, skip_ind);
    end
     
    if is_cell
        obj_update_0{object_id,min(end,layer)} = obj_update;
        obj_illum_0{object_id,min(end,layer)} = obj_illum;
        obj_update = obj_update_0;
        obj_illum = obj_illum_0;
    end
end


function [object_c,object_r] = set_views_gpu_rc(object_c,object_r,proj_c, proj_r,oROI,ind, ind_ok )
% Description: update complex and real object at once 
    import engines.GPU.GPU_wrapper.*

    x =  uint16(oROI{1}(ind,1));
    y =  uint16(oROI{2}(ind,1));

    return_cell = iscell(object_c);

    if ~iscell(object_c); object_c = {Garray(object_c)}; end
    if ~iscell(object_r); object_r = {Garray(object_r)}; end
    if ~iscell(ind_ok); ind_ok = {Garray(uint16(ind_ok))}; end
    proj_c = complex(proj_c); 
    
    set_views_gpu_mex( proj_c, object_c,x,y, ind_ok );
    set_views_gpu_mex( proj_r, object_r,x,y, ind_ok );

    if ~return_cell
       object_c = object_c{1};
       object_r = object_r{1};
    end

end


