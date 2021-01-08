%  FP16.SET simple wrapper around function "halfprecision"
%
% fp_data = set(data)
%
% Inputs: 
%  **data - single precision numerical array
%
% *returns*:
%  ++fp16_data - data are returned as uint16 class which is used as container for fp16 numerical precision 
%
% Recompile command:
% mex -R2018a CXXFLAGS='$CXXFLAGS -fpermissive' 'CFLAGS="\$CFLAGS -fopenmp"' LDFLAGS="\$LDFLAGS -fopenmp"  halfprecision.cpp


%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2018 by Paul Scherrer Institute (http://www.psi.ch)    |
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


function fp16_data = set(data)
    if isscalar(data) && ~isa(data, 'uint16') 
        data = single(data);  % scalar doubles convert to singles directly
    end
    if ~isa(data, 'uint16')
        fp16_data = halfprecision(data); 
    else
        error('Input has to be normal precision')
    end
end
