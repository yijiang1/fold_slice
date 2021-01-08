% GET_VIEWS extract view for each the provided positions and indices 
%
% obj_proj = get_views(object, obj_proj,layer_ids,object_id, indices, cache, scan_ids, skip_ind)
%
% ** object   [Nx_o, Ny_o] array or cells containing object 
% ** obj_proj  [Nx_o, Ny_o, N] preallocated array for the views
% ** layer_ids  id of the solved layer for multilayer ptycho 
% ** object_id  id of the object, ie scan or incoherent mode 
% ** indices   processed positions 
% ** cache     structure with precalculated values to avoid unnecessary overhead
% ** scan_ids   determines to which scan correponds each of the position 
% ** skip_ind   list of indices to be skipped 
%
% returns:
% ++ obj_proj    [Nx_p, Ny_p, N] array with the object views 


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



function obj_proj = get_views(object, obj_proj,layer_ids,object_id, indices, cache, scan_ids, skip_ind)

    import engines.GPU_MS.GPU_wrapper.*
    import engines.GPU_MS.shared.*
    import utils.verbose
    global use_gpu
        
   Np_p = [size(obj_proj,1),size(obj_proj,2)];
    
   if size(obj_proj,3) ~= numel(indices)  || isempty(obj_proj)
        %  in the last iter, number  of items may not be equal to grouping 
        if isempty(obj_proj)
            obj_proj = Gzeros([Np_p, numel(indices)], true);
        else
            % should be tiny bit faster 
            obj_proj = zeros([Np_p, numel(indices)], 'like', obj_proj);
        end
   end
    
   if nargin < 8
       skip_ind = [];
   end
   skip_ind = [skip_ind, cache.skip_ind];
   
   
    if nargin > 6   
        %% !!! call recursivelly -> wrapper for multiscan version !!!!
        
        % get unique IDs of the scans
        if isempty(scan_ids)
            unq_scans = []; 
        elseif all(scan_ids == scan_ids(1))
            unq_scans = scan_ids(1);
        else
            unq_scans = unique(scan_ids);
        end
        if (length(unq_scans)> 1 || length(object) > 1 ) && iscell(object)
           if  ~isempty(use_gpu) && use_gpu && isa(object{object_id(1)}, 'gpuArray') && isa(obj_proj, 'gpuArray') 
               % object_modes > 1 not implemented yet        
               if size(object,1) == 1
                    % shared object or single object 
                    ind_ok{1} = uint16(1:length(indices));
               else
                    for kk = 1:size(object,1)
                        ind_ok{kk} = uint16(find(scan_ids == kk));
                    end
               end
               % feed data directly to the GPU mex without splitting 
               obj_proj = get_views_gpu(object(:,layer_ids),obj_proj,cache.oROI_s{min(end,object_id(1))}, indices, ind_ok);
           else          
               % if GPU not available use this "wrapper" around single
               % set_projection function 
                for kk = unq_scans
                     ind = scan_ids == kk;
                     skip_ind = indices(~ind);  % avoid going through these indices 
                     obj_proj = get_views(object{kk,layer_ids},obj_proj, 1, object_id, indices, cache, [], skip_ind);
                end
           end
           return 
        end
    end
    
    object_id = object_id(1); 
    try
    if iscell(object)
        object = object{object_id, layer_ids};
    end
    catch
        keyboard
    end 
    
    if ~isempty(cache.skip_ind) && ~isempty(skip_ind)
        ind_ok = uint16(find(~ismember(indices, [cache.skip_ind,skip_ind])));  % skip wrong patterns 
    else
        ind_ok = uint16(1:length(indices)); 
    end


    if isempty(obj_proj)
        Np_p = [length(cache.oROI{min(object_id,end)}{1,1}), length(cache.oROI{min(object_id,end)}{1,2})];
        obj_proj = Gzeros([Np_p, length(ind_ok)], true);
    end
    
       
    if isa(object, 'gpuArray')
        %% USE CUDA MEX FOR GPU 
        obj_proj = get_views_gpu(object,obj_proj,cache.oROI_s{min(end,object_id)}, indices, ind_ok);
    else
        %% USE CPU 
        positions = int32([cache.oROI_s{min(end,object_id)}{1}(indices,1), cache.oROI_s{min(end,object_id)}{2}(indices,1)]);
        obj_proj = utils.get_from_3D_projection(obj_proj,object, positions, ind_ok); 
    end
    
end


function obj_proj = get_views_gpu(object, obj_proj, ROI, ind, ind_ok)
    % Description:  simple method to get GPU based projections from object
    import utils.verbose
    % mexcuda  -output +engines/+GPU_MS/get_views_gpu_mex +engines/+GPU/get_views_gpu_mex.cu
    import engines.GPU_MS.GPU_wrapper.*

    
    %% do not return matrices, write directly into obj_proj !!! 
    x =  uint16(ROI{1}(ind,1));
    y =  uint16(ROI{2}(ind,1));

    if ~iscell(object); object = {Garray(complex(object))}; end
    if ~iscell(ind_ok); ind_ok = {Garray(uint16(ind_ok))}; end
    obj_proj = complex(Garray(obj_proj));
    try
        get_views_gpu_mex( obj_proj, object,x,y, ind_ok);
    catch err 
        verbose(0, 'Recompilation of MEX functions ... ')
        if any(strcmp(err.identifier, { 'MATLAB:UndefinedFunction','MATLAB:mex:ErrInvalidMEXFile'}))
            path = replace(mfilename('fullpath'), mfilename, ''); 
            mexcuda('-output', [path,'private/get_views_gpu_mex'], [path, 'private/get_views_gpu_mex.cu'])

            get_views_gpu_mex( obj_proj, object,x,y, ind_ok);
        else
           rethrow(err) 
        end
    end
end

