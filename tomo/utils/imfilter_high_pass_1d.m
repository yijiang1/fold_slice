% IMFILTER_HIGH_PASS_1D applies fft filter along AX dimension that
% removes SIGMA ratio of the low frequencies 
% 
% img = imfilter_high_pass_1d(img, ax, sigma, padding)
%
% Inputs:
%   **img           - ndim filtered image 
%   **ax            - filtering axis 
%   **sigma         - filtering intensity [0-1 range],  sigma <= 0 no filtering 
%   **padding       - pad the array to avoid edge artefacts (in pixels) [default==0]
%   **apply_fft     - if true assume that img is in real space, [default == true]
% *returns*
%   ++img           - highpass filtered image 

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


function img = imfilter_high_pass_1d(img, ax, sigma, padding, apply_fft)

    if nargin < 4, padding = 0; end 
    if nargin < 5, apply_fft = true; end 
    
    Ndims = ndims(img);

    padding = ceil(padding);

    if padding > 0
        pad_vec = zeros(Ndims,1);
        pad_vec(ax) = padding;   % padding only along axis "ax"
        img = padarray(img,pad_vec,'symmetric','both');
    end
    
    Npix = size(img);
    shape = ones(1,Ndims);
    shape(ax) = Npix(ax);
    isReal = isreal(img); 
    
    if apply_fft
        img = fft(img,[],ax);
    end

    x = reshape((-Npix(ax)/2:Npix(ax)/2-1)/Npix(ax), shape);
    
    sigma = 256/(Npix(ax)-2*padding)*sigma; % solution to make the filter resolution independend -> 
                                % -> for different level of scaling the filtered field will look the same 
                                
    if sigma == 0
        % use derivative filter 
        spectral_filter = 2i*pi*(fftshift((0:Npix(ax)-1)/Npix(ax))-0.5);
    else
        spectral_filter = fftshift(exp(1./(-(x.^2)/(sigma)^2)));      
    end
    
    img = bsxfun(@times, img, spectral_filter);

    if apply_fft
        img = ifft(img,[],ax);
    end
    if isReal
        img = real(img);
    end

    if padding > 0
        crop_vec = repmat({':'},Ndims,1);
        crop_vec{ax} = 1+padding:(size(img,ax)-padding);
        img = img(crop_vec{:});
    end

        

end
