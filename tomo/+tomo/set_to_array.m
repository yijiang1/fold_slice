% SET_TO_ARRAY simplified wrapper around CPU-based multithread mex function "add_to_3D_volume_mex"
%
%       set_to_array(full_object, object_block, offset, add_values = false)
%
% Equivalent but faster to matlab command 
%  full_object(:,:, offset + 1:size(object_block,3)) = full_object(:,:,offset + 1:size(object_block,3)) + object_block
%
% Inputs
%   **full_object   - (3D array) array to be added to, note that directly the provided array will be modified, without copying by matlab 
%   **object_block  - (2D/3D array) values to be added to the main structure, can be array or a sharememory class @shm
%   **offset        -(uint scalar) starting position along the 3rd-axis, assuming that full block is loaded
%   **add_values    - (bool) if true, the values will be added otherwise overwritten 
%  *returns*:
%   ++full_object or none - (values are written directly to full_object if MEX access is used)
    
        
%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%*-----------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing language this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS matlab package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.

function full_object = set_to_array(full_object, object_block, offset, add_values)

    if nargin < 4 
        add_values= false; % add values instead of rewritting 
    end
    if isa(object_block,'shm')
        [s,object_block] = object_block.attach(); 
    end
    
    % write back to the full array stored in RAM     
    positions = zeros(size(object_block,3),2); 
    indices = int32(1:size(object_block,3)) + int32(offset(1)); 
    % use a MEX code to speed it up 
    full_object = utils.add_to_3D_projection(object_block, full_object,positions,indices, add_values); 


    if exist('s','var')
        s.free; 
    end
end
    
