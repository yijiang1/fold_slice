% GET_IMG_INT_1 use FFT to integate the image along selected axis -> can be used for phase
% unwrapping
%
%   integer = get_img_int_1D(grad_array, ax)
%
% Inputs:
%   **grad_array   -  phase gradient
%   **ax           -  integration axis 
% *returns*
%   ++integral - 2D stack scalar arrays that has gradients along "ax" close to "grad_array"



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

function integer = get_img_int_1D(grad_array, ax)
    if nargin < 2
        ax = 1;
    end
    
    Np = size(grad_array);    

    if ax == 2
        grad_array = fft(grad_array,[],2);
        xgrid = ifftshift(-fix(Np(2)/2):ceil(Np(2)/2)-1)/Np(2);

        % not sure why, but it seems to need also shift by 1 pixel to make it
        % consistent with gradient 
        X = exp((2i*pi)*xgrid);
        % integration filter
        X = X./ (2i*pi*xgrid);
        X(1) = 0;
        integer = bsxfun(@times, grad_array,X);
        integer = ifft(integer,[],2);

    elseif ax == 1
        grad_array = fft(grad_array,[],1);
        ygrid = ifftshift(-fix(Np(1)/2):ceil(Np(1)/2)-1)/Np(1);

        % not sure why, but it seems to need also shift by 1 pixel to make it
        % consistent with gradient 
        Y = exp((2i*pi)*ygrid');
        % integration filter
        Y = Y./(2i*pi*ygrid');
        Y(1) = 0;
        integer = bsxfun(@times, grad_array,Y);
        integer = ifft(integer,[],1);
    else
        error('Non implemented dimension')
    end
    
        
end
