% INTERPOLATE_LINEAR rescaling based on interp2, faster than utils.interpolateFT,
% works also with GPU 
% Note: for small arrays processed on GPU, utils.interpolateFT can be
% faster due to lower overhead (no for-loop)
%
% img = interpolate_linear(img, scale, method)
%
% Inputs:  
%    **img      -  2D or stack of 2D images 
%    **scale    - scaling factor 
%    **method   - linear (default), cubic, nearest
% *returns*:
%    ++img_out  - 2D or stack of 2D images scaled by factor scale

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



function img_out = interpolate_linear(img, sizeOut, method)

    [Nx, Ny,Nlayers] = size(img);
    
    if all([Nx,Ny] == sizeOut(1:2))
        img_out = img; 
        return
    end
    
    if nargin < 3
        method = 'linear'; 
    end
    
    img_out = zeros([sizeOut(1:2), Nlayers],'like',img);
    for ii = 1:Nlayers
        img_out(:,:,ii) = interp2(img(:,:,ii), linspace(1,Ny,sizeOut(2)),linspace(1,Nx,sizeOut(1))', method);
    end
   
end