%GET_UNIT_LENGTH returns the SI unit for a given length
% [unit,val] = get_length_unit(val)
% Assumes that val is given in [m].
%
%   **val...         length
%
%   returns:
%   ++ unit          SI unit string
%   ++ val           input value converted to SI unit
%
%   EXAMPLE:
%       [unit, val] = get_unit_length(2e-3);
%           unit 
%               'mm'
%           val
%               '2'
%       

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

function [unit,val] = get_unit_length(val)

units = {'m', 'mm', 'um', 'nm', 'pm', 'fm', 'am'};

scl = 0;
while true
    if abs(val)*1e3^(scl)>=1 || scl==length(units)-1
        val = val*1e3^(scl);
        unit = units{scl+1};
        break;
    else
        scl = scl+1;
    end
end


end

