% GARRAY function that helps to wrap the GPU functions for user so that CPU and
% GPU code is identical 
% If gpuDeviceCount > 0 or move_on_GPU == true, the returned array will be 
% moved to GPU, otherwise it will be returned as single 
%
% array = Garray(array, move_on_GPU = true) 
% 
% Inputs: 
%   **array         Ndim array 
%   **move_on_GPU   if true, use GPU if possible 
% retuns:
%   ++array         Ndim array single or gpuArray single 

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


function array = Garray(array, move_on_GPU)
    persistent use_gpu
    if nargin == 2 && ~isempty(move_on_GPU)
        use_gpu = move_on_GPU; 
    elseif isempty(use_gpu)
        use_gpu = gpuDeviceCount > 0;   % always use GPU if not asked otherwise 
    end
    
    if isa(array, 'double') && ~issparse(array)
        %% avoid doubles ... 
        array = single(array);
    end   
   
    if isa(array, 'gpuArray') || ~use_gpu
       return 
    end

    array = gpuArray(array);

end
