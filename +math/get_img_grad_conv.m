% GET_IMG_GRAD_CONV get image gradients along all 3 axis using real space convolution
% 
%  [dX, dY, dZ] = get_img_grad_conv(img, win_size, axis)
%
% Inputs:
%   **img   - stack of images 
%   **win_size - size of the window used for approximation of the FFT gradient 
%   **axis  - direction of the derivative 
% *returns*
%   ++[dX, dY, dZ]  - 3D gradients 

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


function [dX, dY, dZ] = get_img_grad_conv(img, win_size, axis)
    %% get vertical and horizontal gradient of the image 

    if ~isreal(img); error('Not implemented'); end
    ker = get_kernel(win_size);
    if isa(img, 'gpuArray'); ker = gpuArray(ker); end
    if nargin < 3 || any(axis == 2)
        dX = convn(img, reshape(ker,1,[],1), 'same');
    end
    if nargout > 1 || (nargin > 2 && any(axis == 1))
        dY = convn(img, reshape(ker,[],1,1), 'same');
        if nargout == 1; dX = dY; end
    end
    if nargout > 2 || (nargin > 2 && any(axis == 3))
        dZ = convn(img, reshape(ker,1,1,[]), 'same');
        if nargout == 1; dX = dZ; end
    end
end

function ker = get_kernel(win_size)
    N = max(9,2*win_size +1);
    grid = 2i*pi*(fftshift((0:N-1)/(N))-0.5);
    ker = -real(fftshift(fft(grid)))/length(grid); 
    ker = ker( ceil(end/2)+(-ceil(win_size):ceil(win_size)));
    ker = utils.Garray(single(ker));
end