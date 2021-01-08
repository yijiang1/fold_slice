% get vertical and horizontal gradient of the image using FFT
%
%  [dX, dY] = get_img_grad(img, axis, split)
%
% Inputs:
%   **img   - stack of images 
%   **split - split for GPU fft_partial 
%   **axis  - direction of the derivative 
% *returns*
%  ++[dX, dY]  - image gradients 

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
% 



function [dX, dY] = get_img_grad(img, axis, split)
    import math.* 
    
    if nargin < 3 || isempty(split)
        split = 1;
    end
    isReal = isreal(img); 
    Np = size(img);
    if nargin < 2 || any(axis == 2)
        X = 2i*pi*ifftshift(-fix(Np(2)/2):ceil(Np(2)/2)-1)/Np(2);
        dX = bsxfun(@times,fft_partial(img,2,1,split,false),X);
        dX = fft_partial(dX,2,1,split,true);
        if isReal; dX = real(dX);end 
    end
    if nargout == 2 || (nargin > 1 && any(axis == 1))
        Y = 2i*pi*ifftshift(-fix(Np(1)/2):ceil(Np(1)/2)-1)/Np(1);
        dY = bsxfun(@times, fft_partial(img,1,2,split,false),Y.');
        dY = fft_partial(dY,1,2,split,true);
        if isReal; dY = real(dY);end 
        if nargout == 1; dX = dY; end
    end
end

