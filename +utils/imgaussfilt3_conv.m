% IMGAUSSFILT3_CONV apply gaussian smoothing along all three dimensions using convolution, 
% faster than matlab alternative
%
%
%  A = imgaussfilt3_conv(A,sigma)
%
% Inputs: 
%   **A           3D volume to be smoothed 
%   **sigma       gaussian smoothing constant, scalar or use vector for anizotropic kernel smoothing
% returns: 
%   ++A           smoothed volume 

        

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



function X = imgaussfilt3_conv(X, filter_size)
    %% faster equivalent to the imgaussfilt3 in matlab 
        
    shape_0 = {[], 1,1};
    for ax = 1:3
        if filter_size(min(end,ax))  == 0 
            continue
        end
        if ax == 1 || filter_size(min(end,ax-1)) ~= filter_size(min(end,ax)) 
            ker = get_kernel(filter_size(min(end,ax)) , class(X));
        end
        shape = circshift(shape_0, ax-1);       
        X = convn(X, reshape(ker,shape{:}), 'same');       
    end
    
end

function ker = get_kernel(filter_size, class)
    
    grid = (-ceil(2*filter_size):ceil(2*filter_size)) / filter_size;
    ker = exp(-grid.^2);
    ker = ker / sum(ker); 
    if isa(class, 'gpuArray')
        ker = gpuArray(single(ker));
    end
end
