% FIND_SHIFT_FAST_1D uses cross-correlation to find shift between 1D patterns o1 and
% o2, if the patterns are 2D, perform the search along the axis `ax`
%
% shift = find_shift_fast_1D(o1, o2,  ax, sigma)
%
% Inputs:
%   **o1            - aligned array 1D/2D (will be aligned along 1st axis)
%   **o2            - template for alignment 1D or 2D 
%   **ax            - perform search along this axis 
%  *optional*
%   **sigma         - filtering intensity [0-1 range],  sigma <= 0 no filtering, recommended sigma < 0.05 
%   **padding       - pading [in pixels] the provided array by zeros, prevent circular boundary condition in FFT, default = 0 
% *returns*
%   ++shift         - (vector) displacement of the 1D/2D arrays  

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



function shift = find_shift_fast_1D(o1, o2, ax, sigma, padding)

    if nargin < 3
        ax = 2; 
    end
    if nargin < 4
        sigma = 0;
    end
    if nargin < 5
        padding = 0; 
    else
        padding = ceil(padding/2)*2;
    end
    
    max_shift = size(o1,ax)/3;  % avoid too large corrections !!!
    
    Ndims = ndims(o1);
   
    if ax ~= 1 
        error('FIXME: Not tested axis')
    end
    
    %% symmetrize before spectral filtering !!
    o1 = cat(ax, o1, flipud(o1)); 
    o2 = cat(ax, o2, flipud(o2)); 
    Npix = size(o1);

    shape = ones(1,Ndims);
    shape(ax) = Npix(ax);

    if sigma > 0
        %% high pass filter 
        o1 = fft(o1, [],ax);
        o2 = fft(o2, [],ax);
        x = reshape((-Npix(ax)/2+1:Npix(ax)/2)/Npix(ax), shape);
        spectral_filter = fftshift(exp(1./(-(x.^2)/(sigma^2))));
        spectral_filter(floor(end/2+[-3:3])) = 0; %% remove some strange artefacts 
        o1 = bsxfun(@times, o1, spectral_filter);
        o2 = bsxfun(@times, o2, spectral_filter);
        o1 = ifft(o1, [],ax);
        o2 = ifft(o2, [],ax);
    
    end

    % remove symetrization  !!
    o1 = o1(1:end/2,:);
    o2 = o2(1:end/2,:);
    o1 = padarray(o1, padding/2, 'both');
    o2 = padarray(o2, padding/2, 'both');
    
    
    Npix = size(o1);
    shape = ones(1,Ndims);
    shape(ax) = Npix(ax);

    %% remove edge issues (after symetrized filtering )
    spatial_filter = reshape(tukeywin(prod(shape)), shape);
    o1 = bsxfun(@times, o1, spatial_filter);
    o2 = bsxfun(@times, o2, spatial_filter);
    
    o1 = fft(o1, [],ax);
    o2 = fft(o2, [],ax);

    %% cross-correlation 
    xcorrmat = abs(ifft(o1.*conj(o2),[],ax)); 
    %% 1D fftshift 
    xcorrmat = circshift(xcorrmat, floor(Npix(ax)/2), ax); 
    if ax == 2; error('FIXME: Not tested axis'); end 
    
    % choose only optimim withing reduced range 
    xcorrmat([1:ceil(end/2-max_shift), ceil(end/2+max_shift):end],:) = 0; 
    
    %% take only small region around maximum 
    WIN = 10; 
    kernel_size = [1,1];
    kernel_size(ax) = WIN; 
    mask = conv2(single(bsxfun(@eq, xcorrmat, max(xcorrmat,[],ax))), ones(kernel_size), 'same'); 
    xcorrmat(~mask) = nan; 
    xcorrmat = max(0, bsxfun(@minus, xcorrmat, min(xcorrmat,[],ax)));
    xcorrmat(~mask) = 0; 
    xcorrmat = bsxfun(@times, xcorrmat, 1./max(xcorrmat,[],ax)).^4;
    
    
    %% find center of mass 
    MASS = sum(xcorrmat,ax);
    grid = reshape(1:Npix(ax),shape); 
    shift = sum(bsxfun(@times, xcorrmat, grid),ax) ./ MASS - floor(Npix(ax)/2)-1;
    
end
