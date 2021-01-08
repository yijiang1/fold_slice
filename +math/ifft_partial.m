% IFFT_PARTIAL apply fft only on smaller blocks (important for GPU)
% 
%  x = ifft_partial(x,fft_axis,split_axis, split)
% 
% Inputs: 
%   **x         - input array 
%   **fft_axis  - axis along which is performed FFT 
%   **split_axis - axis along which is the array split 
%
%  *optional*
%   **split - number of blocks to split the array before FFT to save the memory 
% 
% *returns*
%  ++x       ifft transformed array 

    
    
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


function x = ifft_partial(x,fft_axis,split_axis, split)

    if nargin < 4
        split = []; 
    end
 
    x = math.fft_partial(x,fft_axis,split_axis, split, true); 

    
end