% IMROTATE_AX_FFT fft-based image rotation for a stack of images along given axis 
% based on "Fast Fourier method for the accurate rotation of sampled images", Optic Communications, 1997
% img_stack = imrotate_ax_fft(img, theta, ax)
% 
% Inputs: 
%   **img       - stacked array of images to be rotated 
%   **theta     - rotation angle 
% *optional*
%   **axis      - rotation axis (default=3)
% returns:
%   ++img       - rotated image 

    
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


function img = imrotate_ax_fft(img, theta, axis)
    if all(theta == 0) || isempty(img); return ; end

    if nargin  < 3
        axis = 3;
    end
    
    isReal = isreal(img); 
    
    if axis == 1
       img = permute(img, [3,2,1]); theta = - theta;
    elseif axis == 2
       img = permute(img, [1,3,2]);
    end
    
    angle_90_offset = round(theta/90);

    if angle_90_offset ~= 0
        img = rot90(img, angle_90_offset); 
        theta = theta - 90*angle_90_offset;
    end
    
    if theta == 0; return ; end


    [M, N, ~] = size(img); 

    % make possible to rotate each slice with different angle 
    theta = reshape(theta,1,1,[]) * ones(1,'like',img); % move to GPU if needed 
    xgrid = (ifftshift(-fix(M/2):ceil(M/2)-1)'/M); 
    ygrid = (ifftshift(-fix(N/2):ceil(N/2)-1) /N); 
    Mgrid = (1:M)'-floor(M/2)-0.5; % the 0.5px offset is important to make the rotation equivalent to matlab imrotate
    Ngrid = (1:N) -floor(N/2)-0.5; 
    
    if isa(theta, 'gpuArray')
        [M1, M2] = arrayfun(@aux_fun, theta, xgrid, ygrid, Mgrid, Ngrid); 
    else
        [M1, M2] = aux_fun(theta, xgrid, ygrid, Mgrid, Ngrid); 
    end 
   
    
    % rotate images by a combination of shears 
    img=ifft(fft(img,[],2).*M1,[],2);
    img=ifft(fft(img,[],1).*M2,[],1);
    img=ifft(fft(img,[],2).*M1,[],2);

    if isReal
        img = real(img); 
    end

    if axis == 1
       img = permute(img, [3,2,1]);
    elseif axis == 2
       img = permute(img, [1,3,2]);
    end

  
end

% auxiliarly function to be used for GPU kernel merging
function [M1, M2] = aux_fun(theta, xgrid, ygrid, Mgrid, Ngrid)

    
    % based on "Fast Fourier method for the accurate rotation of sampled images", Optic Communications, 1997
    Nx = -sind(theta) .* xgrid; 
    Ny = tand(theta/2).* ygrid;
    
    M1 = exp(-2i*pi*Mgrid.*Ny); 
    M2 = exp(-2i*pi*Ngrid.*Nx); 

end