%  GET_FROM_ARRAY simplified wrapper around CPU-based multithread mex function "get_from_3D_projection"
%
%  object_block = get_from_array(full_object, object_block, indices, offset)
%
%
% Inputs:  
%    **full_object          - (3D array) large volume block where data will be loded from 
%    **object_block         - (3D array) smaller volume block into which the selected data from full_object will be added 
%    **indices              - list of slices along 3rd axis that are written to object_block, starting from 0
%    **offset               - (numel(indices)x2 nonnegative integers) offset from (1,1) corner
% *returns*
%    ++object_block         - loaded small block of the full_object 
% Example: 
%   function is equivalent (but much faster) to 
%   for ii = 1:size(object_block,3)
%       object_block(:,:,ii) = full_object(offset(ii,1)+1:size(object_block,1) , offset(ii,2)+1:size(object_block,2),ii);
%   end
    
    
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

function object_block = get_from_array(full_object, object_block, indices, offset)
    if nargin < 4
        offset = [0,0];
    end
    
    Nind = length(indices); 
    if isempty(object_block)
        object_block = zeros(size(full_object,1), size(full_object,2), Nind, 'like',full_object); 
    end
    if isa(object_block,'shm')
         object_block.allocate(zeros(size(full_object,1), size(full_object,2), Nind, 'like',full_object)); 
        [s,object_block] = object_block.attach();
    end
    
    positions = zeros(size(object_block,3),2,'int32') + int32(offset);
    
    % mex custom made fast get function 
    utils.get_from_3D_projection(object_block,full_object,positions, int32(indices));
    
    
    if exist('s','var')
        s.detach; 
        object_block = s; 
    end
end
