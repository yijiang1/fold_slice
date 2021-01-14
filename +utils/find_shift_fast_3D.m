% FIND_SHIFT_FAST_3D uses crosscorelation to find shift between o1 nd
% o2 patterns in 3D space 
%
% shift = find_shift_fast_3D(o1, o2, sigma, apply_fft)
%
% Inputs:
%   **o1            - aligned array  3D - return only single shift vector [x,y,z]
%   **o2            - template for alignement  3D
%   **sigma         - filtering intensity [0-1 range],  sigma <= 0 no filtering, recommended sigma < 0.05 
%   **apply_fft     - if false, assume o1 and o2 to be already in fourier domain 
% *returns*
%   ++shift         - displacement of the 3D volumes 

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


function shift = find_shift_fast_3D(o1, o2, sigma, apply_fft)   
    

    import math.*
    
    assert(ndims(o1) == 3, 'Inputs has to be 3D matrix')
    assert(ndims(o2) == 3, 'Inputs has to be 3D matrix')
    
   if nargin < 4 
       apply_fft = true;
   end
   if nargin < 3
       sigma = 0.01;
   end
  
   
   if apply_fft

        [nx, ny, nz] = size(o1);

        % suppress edge effects of the registration procedure
        spatial_filter = tukeywin(nx,0.5) * tukeywin(ny,0.5)' .* reshape(tukeywin(nz,0.5),1,1,[]); 
        
        o1 = bsxfun(@times, o1, spatial_filter);
        o2 = bsxfun(@times, o2, spatial_filter);

        clear spatial_filter
        
        o1 = fftn(o1);
        o2 = fftn(o2);
        
   end
   
    [nx, ny, ~] = size(o1);

    if sigma > 0 
        % remove low frequencies 
        [X,Y,Z] = meshgrid( (-nx/2:nx/2-1)/nx, (-ny/2:ny/2-1)/ny, (-nz/2:nz/2-1)/nz);
        spectral_filter = fftshift(exp(1./(-(X.^2+Y.^2+Z.^2)/sigma^2)));
        o1 = bsxfun(@times, o1, spectral_filter);
        o2 = bsxfun(@times, o2, spectral_filter);  
        clear spectral_filter
    end
  
    
    % fast subpixel cross correlation 
    xcorrmat = fftshift(abs(ifftn(o1.*conj(o2)))); 
    

    %% take only small region around maximum 

    WIN = 5; 
    kernel_size = [WIN,WIN,WIN];
    xcorrmat = xcorrmat / max(xcorrmat(:)); 
    xcorrmat = xcorrmat .* convn(xcorrmat == 1, ones(kernel_size,'single'), 'same'); 
    [x,y,z] = find_center_fast(xcorrmat.^2);
    shift = [x,y,z];

end

function [x,y,z] = find_center_fast(xcorrmat)
    MASS = squeeze(sum(xcorrmat(:)));
    [N,M,O] = size(xcorrmat); 
    x = squeeze(sum(sum(sum(xcorrmat .* (1:M),1)))) ./ MASS - floor(M/2)-1;
    y = squeeze(sum(sum(sum(xcorrmat .* (1:N)',2)))) ./ MASS - floor(N/2)-1;
    z = squeeze(sum(sum(sum(xcorrmat .* reshape(1:O,1,1,[]),3) ))) ./ MASS - floor(O/2)-1;
end
