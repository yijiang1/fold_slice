% BINNING_3D - bin data along first three axis, 
%  x = binning_3D(x, binning, centered)
%
%  Inputs: 
%   **x - original array, upsampling will be performed only along the first three axis, array size has to be dividable by binning size 
%   **binning - scalar or (3,1) array,  positive integer binning factor 
%  *optional*
%   **centered - default false, shift the binning by binning/2 offset 
%
%   returns: 
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



function x = binning_3D(x, binning, centered)
    if all(binning == 1); return ; end 
    if nargin < 3
        centered = false; 
    end
    if ismatrix(x)
        x = utils.binning_2D(x,binning, centered);
        return
    end
    binning = reshape(binning,1,[]);

    assert(ndims(x) == 3, 'Input has to be 3D array')
    assert(all(mod(size(x),binning)==0), 'Array has to splitable by binning')
    assert(all(size(x) >binning ), 'Array has to larger than binning')

    Npix = size(x); 
    if isscalar(binning)
        binning = repmat(binning,3,1); 
    end
    
    if centered
        % make the bins centered 
        x = x(ceil(binning(1)/2):end-ceil(binning(1)/2)-1, ceil(binning(2)/2):end-ceil(binning(2)/2)-1,ceil(binning(3)/2):end-ceil(binning(3)/2)-1); 
        Npix(1:3) = Npix(1:3) - reshape(binning,1,[]); 
    end
        
    x = reshape(x,binning(1), Npix(1)/binning(1), binning(2), Npix(2)/binning(2), binning(3), Npix(3)/binning(3) ); 
    x = squeeze(sum(sum(sum(x,1),3),5)); 
    x = x / prod(binning); 

end