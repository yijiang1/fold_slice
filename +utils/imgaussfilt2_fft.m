% IMGAUSSFILT2_FFT apply gaussian smoothing along all three dimensions
% faster than matlab version 
%
%  A = imgaussfilt2_fft(A,sigma)
%
% Inputs: 
%   **A   3D volume to be smoothed 
%   **sigma   gaussian smoothing constant 
%   **split   3x1 int vector to split the volume and save memory 
% *returns*: 
%   ++A    smoothed volume 

        

    
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



function A = imgaussfilt2_fft(A,sigma, split)
    % gaussian blurring along first 2 dimension 
    import math.*
    if isscalar(sigma) && sigma == 0
        return
    end
    if nargin < 3
        split = 1; 
    end
    
    isReal = isreal(A); 
    
    Npx = size(A);
    
    A = fft2_partial(A, split);

    for dim=1:2
        if isscalar(sigma)
            grid = single((-Npx(dim)/2:Npx(dim)/2-1));
            ker = exp(-grid.^2/ sigma^2)';
        elseif isvector(sigma)
            % use user given kernel , assume splitable 1D kernel 
            ker = zeros(Npx(dim),1,'like', A); 
            Ns = length(sigma);
            ker(ceil(Npx(dim)/2)+[-floor(Ns/2):ceil(Ns/2)-1]) = sigma;
        else
            error('N-dim kernel not implemented')
        end
        ker = ker / sum(ker);
        ker = fft(ker,[],1);

        ker_shape = ones(1,2);
        ker_shape(dim) = Npx(dim);
        B = reshape(ker,ker_shape);
        if isa(A, 'gpuArray'); B = gpuArray(B); end
        A = A.*B;
    end
    clear B
    
    A = ifft2_partial(A, split);
    
    if isReal
        A = real(A);
    end
    
    % Im not sure why, but the output needs to be fftshifted
    A = fftshift_2D(A);
        
end