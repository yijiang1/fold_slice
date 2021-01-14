% INTERPOLATEFT_CENTERED Perform FT interpolation of provided stack of images using FFT so that
% the center of mass is not modified after the resolution change
% This function is critical for subpixel accurate up/down sampling
%
% imout = interpolateFT_centered(im,downsample,interp_sign)
%
% Inputs:
%   **im                - Input complex  2D array or stacked 3D array 
%   **Npix_new          - (2x1 vector) Size of the interpolated array
%   **interp_sign       - +1 or -1, sign that adds extra 1px shift. +1 is needed
%                           if interpolation is used to downsample phase gradient which is used for
%                           unwrapping, otherwise use -1
% *returns*:
%   ++imout   - Output complex image
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

function [ img ] = interpolateFT_centered(img,Np_new, interp_sign)
    import utils.*
    import math.*

    Np = size(img); 
    Np_new = 2+Np_new; 
    isReal = isreal(img); 

    scale = prod((Np_new-2)) / prod(Np(1:2));
    downsample = ceil(sqrt(1/scale)); 

    if isa(img, 'gpuArray')
        scale = Garray(scale); 
    end
    % apply the padding to account for boundary issues 
    img = padarray(img, double([downsample,downsample]), 'symmetric' ,'both'); 

    % go to the fourier space 
    img = fft2(img);

    % apply +/-0.5 px shift
    img = imshift_fft(img, interp_sign*-0.5, interp_sign*-0.5, false);

    % crop in the Fourier space (can be speeded up similarly to example in utils.interpolateFT_ax )
    img = ifftshift_2D(crop_pad(fftshift_2D(img), Np_new));

    % apply -/+0.5 px shift in the cropped space 
    img = imshift_fft(img,  interp_sign*0.5, interp_sign*0.5, false);

    % return to the real space 
    img = ifft2(img);

    % scale to keep the average constant 
    img = img*scale;

    % remove the padding
    img = img(2:end-1, 2:end-1,:); 

    if isReal
        img = real(img); % preserve complexity 
    end
end
