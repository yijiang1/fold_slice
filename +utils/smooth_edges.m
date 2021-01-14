%   SMOOTH_EDGES takes stack of 2D images and smooths boundaries to avoid sharp edge artefacts during imshift_fft 
% 
%   img = smooth_edges(img, win_size, dims)
%
%   Inputs:
%       **img - 2D stacked array, smoothing is done along first two dimensions 
%       **win_size - size of the smoothing region, default is 3 
%       **dims - list of dimensions along which will by smoothing done 
%   Outputs: 
%       ++img - smoothed array 
    
% *-----------------------------------------------------------------------*
% |                                                                       |
% |  Except where otherwise noted, this work is licensed under a          |
% |  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
% |  International (CC BY-NC-SA 4.0) license.                             |
% |                                                                       |
% |  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
% |                                                                       |
% |       Author: CXS group, PSI                                          |
% *-----------------------------------------------------------------------*
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
% 
% 

function img = smooth_edges(img, win_size, dims)

    if nargin < 3
        dims = [1,2];  % default is smooth along first 2 dimensions 
    end
    
    if nargin < 2
        win_size = 5; 
    end
   
    try
        Npix = size(img);
        for i = dims
            ind = {':',':',':'};
            if Npix(i) <= 2*win_size ; continue; end
            win_size = max(win_size, 3); 
            % get  indices of the edge regions 
            ind{i} = [Npix(i)-win_size+1:Npix(i),1:win_size];
            ker_size = [1,1]; 
            ker_size(i) = win_size;
            img_tmp = img(ind{:}); 
            kernel = reshape(gausswin(win_size,2.5), ker_size); 
            % smooth across the image edges 
            img_tmp = convn(img_tmp, kernel, 'same');
            % avoid boundary issues from convolution 
            boundary_shape = [1,1]; 
            boundary_shape(i) = length(ind{i}); 
            img_tmp = bsxfun(@rdivide, img_tmp, conv(ones(boundary_shape), kernel, 'same'));
            img(ind{:}) = img_tmp; 
        end
    catch err
        warning('Smooth edges failed: %s', err.message)
    end

end

function w = gausswin(L, a)
    % Compute window according to [1]
    N = L-1;
    n = (0:N)'-N/2;
    w = exp(-(1/2)*(a*n/(N/2)).^2);
end

