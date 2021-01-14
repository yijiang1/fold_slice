% BINNING_2D - bin data along first two axis, 
%  x = binning_2D(x, binning, centered)
%
%  Inputs: 
%   **x - original array, upsampling will be performed only along the first two axis, array size has to be dividable by binning size 
%   **binning - scalar or (2,1) array,  positive integer binning factor 
%  *Optional*:
%   **centered - default false, shift the binning by binning/2 offset 
%
%   Outputs:
%   ++x - binned array 

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




function x = binning_2D(x, binning, centered)
    
    if all(binning <= 1); return ; end 

    if nargin < 3 
        centered = false; 
    end
    Npix = [size(x,1), size(x,2), size(x,3)]; 
    if isscalar(binning); binning = repmat(binning, 1,2); end
    binning = reshape(binning,1,[]);
    
    if all(Npix(1:2) >= binning(:))
        % faster but less general version 
        if centered
            % it will be slower due to memory copy 
            x = x(ceil(binning(1)/2):end-ceil(binning(1)/2)-1, ceil(binning(2)/2):end-ceil(binning(2)/2)-1,:); 
            Npix(1:2) = Npix(1:2) - binning; 
        end
        if any(~math.isint(Npix(1:2)./binning))
           % is the array cannot be easily split for binning, crop it 
           % it will be slower due to memory copy 
           Npix(1:2) = floor(Npix(1:2)./binning) .* binning;
           x = x(1:Npix(1),1:Npix(2),:); 
        end
        x = reshape(x,binning(1), Npix(1)/binning(1), binning(2), Npix(2)/binning(2), Npix(3)); 
        x = squeeze(sum(sum(x,1),3)); 
        x = x / prod(binning(1:2)); 
        if centered
           x = padarray(x, [1,1], 'replicate', 'post');  % account for the removed pixels to keep the size
        end
    else
        x = convn(single(x), ones(binning, 'single'), 'same'); 
        norm = binning.^2; 
        ind = {ceil(binning(1)/2):binning(1):Npix(1), ceil(binning(2)/2):binning(2):Npix(2)};
        % avoid issues with void dimensions 
        for i = find(Npix == 1)
            ind{i} = ':'; 
            norm = norm / binning(i);  %% avoid summing up by convolution
        end
        x = x(ind{:},:) / norm; 
    end
end