% IMSHEAR_FFT fft-based image shearing function for a stack of images along given axis 
% 
% img_stack = imshear_fft(img_stack, theta, shear_axis)
%
% Inputs: 
%   **img_stack  - stack of 2D images 
%   **theta      - shear angle, scalar 
%   **shear_axis - image axis along which the image will be shared 
% *returns*: 
%   ++img        - shreared image 
%


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


function img = imshear_fft(img, theta, shear_axis)
    
    if theta == 0; return ; end
    
    assert(any(shear_axis == [1,2]), 'Shear axis has to be 1 or 2' )
    
    isReal = isreal(img); 
    
    if abs(theta) > 45
        error('Out of valid angle range [-45,45], use rot90 to get into the valid range')
    end
    
    [M, N, ~] = size(img); 
    theta = reshape(theta,1,1,[]); % allow different theta for each slice 

    Nx = -sind(theta) .* ifftshift(-fix(M/2):ceil(M/2)-1)/M; 
    Ny = tand(theta/2).* ifftshift(-fix(N/2):ceil(N/2)-1)/N;
    Mgrid = 2i*pi*((1:M)'-floor(M/2)) * ones(1,'like',img); 
    Ngrid = 2i*pi*((1:N)'-floor(N/2)) * ones(1,'like',img); 
    
    % rotate images by a combination of shears 
    switch shear_axis
        case 1, img=ifft(fft(img,[],2).*exp(-Mgrid.*Ny), [],2);
        case 2, img=ifft(fft(img,[],1).*exp( Ngrid.*Nx)',[],1);
        otherwise
            error('Shear axis has to be 1 or 2')
    end

    if isReal
        img = real(img); 
    end
   
end