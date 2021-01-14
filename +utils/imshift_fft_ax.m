% IMSHIFT_FFT_AX  will apply subpixel shift that can be different for each 
% frame along one dimension only 
% If apply_fft == false, then images will be assumed to be in fourier space 
%
% Inputs:
%   **img - inputs ndim array to be shifted along ax-th dimension
%   **ax  - axis along which the array will be shifted 
%   **shift - Nx1 vector of shifts, positive direction is up
%   **apply_fft = false - if the img is already after fft, default is false
% *returns*: 
%   ++img - shifted image  / volume 


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


function img = imshift_fft_ax(img, shift, ax, apply_fft)


    if nargin < 4
        apply_fft = true;
    end
    if all(shift == 0)
        return
    end
    
    isReal = isreal(img); 
    
    Npix = size(img);
    
    if ndims(img) == 3
        Np = [1,1,Npix(3)];
    else
        Np = Npix;
        Np(ax) = 1; 
    end

    Ng = ones(1,3);
    if ax > ndims(img)
        Npix(ax) = 1; 
    end
    
    Ng(ax) = Npix(ax);

    
    
    
    if isscalar(shift)
         shift = shift .* ones(Np);
    end
    
    grid = ifftshift(-fix(Npix(ax)/2):ceil(Npix(ax)/2)-1)/Npix(ax);
      
    X = bsxfun(@times, reshape(shift,Np), reshape(grid,Ng));
    X =  exp((-2i*pi)*X);
    
    
    if apply_fft
        img = math.fft_partial(img, ax, 1+mod(ax, ndims(img)) );
    end
        
    img = bsxfun(@times, img,X);
   
    if apply_fft
        img = math.ifft_partial(img, ax, 1+mod(ax, ndims(img)) );
    end
        
    if isReal
        img = real(img);
    end

end

