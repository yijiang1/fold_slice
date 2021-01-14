%   FIND_IMG_ROTATION_2D  find object rotation that provides in projection most sparse features 
%
% [angle] = find_img_rotation_2D(img) 
%
% Inputs: 
%  **img - 2D image to be rotated 
% *returns*: 
%   ++angle - optimal rotation angle in degrees 

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


function  [angle_fine] = find_img_rotation_2D(img, max_range)
    import math.argmin
    if nargin < 2
        max_range = [-22.5,22.5];
    end
    
    %  grid search first => avoid local minimums
    
    test_img = abs(img);
    test_img = max(0, test_img - median(test_img(:)));
    N = 50;
    score  = zeros(N,1);
    alpha_range = linspace(max_range(1),max_range(end), N);
    for i = 1:N
        score(i) = gather(get_score(test_img, alpha_range(i)));
    end
    
    alpha_range = alpha_range(argmin(score)) + (-1:0.1:1);
    clear score
    for i = 1:length(alpha_range)
        score(i) = gather(get_score(test_img, alpha_range(i)));
    end
    angle = alpha_range(argmin(score));
    
    angle_fine = fminsearch(@(x)get_score(test_img, x), angle, struct('TolX', 1e-4)); 
    

    if isa(angle, 'gpuArray')
        angle = gather(angle); 
    end
    fprintf('Optimal image rotation: %.3g°\n', angle_fine)
  
end

function score = get_score(data, angle)
    Npix = size(data);
    [X,Y] = meshgrid(-ceil(Npix(2)/2):floor(Npix(2)/2)-1,-ceil(Npix(1)/2):floor(Npix(1)/2)-1);
    data = data .* (X.^2 / (Npix(2)/2)^2 +Y.^2/(Npix(1)/2)^2 < 1/2); 
    data = data - utils.imgaussfilt2_fft(data,5); 
    
    data = utils.imrotate_ax_fft(data, angle, 3);

    data = data(ceil(end*0.1):floor(end*0.9), ceil(end*0.1):floor(end*0.9));   
    data = (abs(math.fftshift_2D(fft2(data))));
    score =  -mean([sparseness(nanmean(data,1)), ...
        sparseness(nanmean(data,2))]);
    score = gather(score);
end

function spars = sparseness(x)
    %Hoyer's measure of sparsity for a vector
    % from  scipy.linalg import norm

    order_1 = 1;
    order_2 = 2;
    x = x(:);
    sqrt_n = sqrt(length(x));
    spars = (sqrt_n - norm(x, order_1) / norm(x, order_2)) / (sqrt_n - order_1);
end