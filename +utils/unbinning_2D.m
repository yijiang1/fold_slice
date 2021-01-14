% UNBINNING_2D - perform reverse opetation to binning, is replication if
% elements for given upsampling factor. It is very fast method for data
% upsampling
%  x = unbinning_2D(x, upsample)
%  Inputs: 
%   **x - original array, upsampling will be performed only along the first two axis
%   **upsample - scalar or (2,1) array,  positive integer upsampling factor 
%
%   returns:  
%   ++x - upsampled array 


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




function x = unbinning_2D(x, upsample)
    
    if all(upsample <= 1); return ; end 
    
    Npix = [size(x,1), size(x,2), size(x,3)]; 
    if isscalar(upsample); upsample = repmat(upsample, 1,2); end
    upsample = reshape(upsample,1,[]);
    
    x = reshape(x, [1,Npix(1),1,Npix(2),Npix(3)]); 
    x = repmat(x, [upsample(1),1,upsample(2),1,1]); 
    x = reshape(x,Npix(1)*upsample(1), Npix(2)*upsample(2), Npix(3)); 
    
end