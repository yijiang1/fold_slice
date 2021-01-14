% FIND_SHIFT_FAST_2D uses crosscorelation to find shift between o1 nd
% o2 patterns in 3D space 
%
% shift = find_shift_fast_2D(o1, o2, sigma, apply_fft)
%
% Inputs:
%   **o1            - aligned array 2D or 3D, (for stack of images, alignment is done along 3rd axis)
%   **o2            - template for alignment 2D or 3D
%   **sigma         - filtering intensity [0-1 range],  sigma <= 0 no filtering, recommended sigma < 0.05 
%   **apply_fft     - if false, assume o1 and o2 to be already in fourier domain 
% *returns*
%   ++shift         - displacement of the 2D volumes 

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


function shift = find_shift_fast_2D(o1, o2, sigma, apply_fft, method)   
    

    import math.*
    
    
   if nargin < 4 
       apply_fft = true;
   end
   if nargin < 3
       sigma = 0.01;
   end
   if nargin < 5
       method = 'full_range';
   end
   
   if apply_fft

        [nx, ny, ~] = size(o1);

        % suppress edge effects of the registration procedure
        spatial_filter = tukeywin(nx,0.5) * tukeywin(ny,0.5)'; 
        
        o1 = bsxfun(@times, o1, spatial_filter);
        o2 = bsxfun(@times, o2, spatial_filter);

        o1 = fft2(o1);
        o2 = fft2(o2);
   end
   
    [nx, ny, ~] = size(o1);

    if sigma > 0 
        % remove low frequencies 
        [X,Y] = meshgrid( (-nx/2:nx/2-1)/nx, (-ny/2:ny/2-1)/ny);
        spectral_filter = fftshift(exp(1./(-(X.^2+Y.^2)/sigma^2)))';
        o1 = bsxfun(@times, o1, spectral_filter);
        o2 = bsxfun(@times, o2, spectral_filter);       
    end
  
    
    % fast subpixel cross correlation 
    xcorrmat = fftshift_2D(abs(ifft2(o1.*conj(o2)))); 
    

    % %% just for testing 
    % subplot(3,1,1)
    % imagesc(abs(fft2(o1(:,:,1)))); axis off image; colormap bone 
    % subplot(3,1,2)
    % imagesc(abs(fft2(o2(:,:,1)))); axis off image; colormap bone 
    % subplot(3,1,3)
    % imagesc(xcorrmat(:,:,1)); axis off image; colormap bone 
    % drawnow 
    % pause(0.1)

    switch method
        case 'full_range'
            %% take only small region around maximum 
            WIN = 5; 
            kernel_size = [WIN,WIN];
            % convolution may be quite slow ?
            mask = convn(single(bsxfun(@eq, xcorrmat, max2(xcorrmat))), ones(kernel_size,'single'), 'same'); 
            xcorrmat(~mask) = nan; 
            xcorrmat = max(0, bsxfun(@minus, xcorrmat, min2(xcorrmat)));
            xcorrmat(~mask) = 0; 
            xcorrmat = bsxfun(@times, xcorrmat, 1./max2(xcorrmat)).^2;

            %% get CoM of  the central  peak only !!, assume a single peak 
            xcorrmat = max(0, xcorrmat - 0.5).^2;
            [x,y] = find_center_fast(xcorrmat);
            shift = [x,y];

        case 'limited_range'
            % second option: assume that the shifts are only small, it is faster
            mxcorr = mean(xcorrmat,3); 
            [m,n] = find(mxcorr == max(mxcorr(:)));    

            MAX_SHIFT = 10;  % +-10px search 
            MAX_SHIFT_X = min(floor(nx/2-0.5),MAX_SHIFT); 
            MAX_SHIFT_Y = min(floor(ny/2-0.5),MAX_SHIFT); 

            xrange = (-MAX_SHIFT_X:MAX_SHIFT_X); 
            yrange = (-MAX_SHIFT_Y:MAX_SHIFT_Y); 

            idx = { m + xrange,n+yrange,':'};
            xcorrmat = xcorrmat(idx{:});
            MAX = max(max(xcorrmat));

            xcorrmat = bsxfun(@times, xcorrmat, 1. / MAX);

            %% get CoM of  the central  peak only !!, assume a single peak 
            xcorrmat = max(0, xcorrmat - 0.5).^2;
            [x,y] = find_center_fast(xcorrmat);
            shift = [x,y]+[n,m]-floor([ny,nx]/2)-1;
    end
    
    if any(isnan(gather(shift)))
        keyboard
    end
        
end

function [x,y,MASS] = find_center_fast(xcorrmat)
    MASS = squeeze(sum(sum(xcorrmat)));
    [N,M,~] = size(xcorrmat); 
    x = squeeze(sum( bsxfun(@times, sum(xcorrmat,1), 1:M), 2)) ./ MASS - floor(M/2)-1;
    y = squeeze(sum(bsxfun(@times, sum(xcorrmat,2), (1:N)'),1)) ./ MASS - floor(N/2)-1;
end
