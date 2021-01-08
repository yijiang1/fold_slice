%DOUBLE2INT Convert structure values from double to int if precision can be
%preserved.
% 
% EXAMPLE: 
%   p.value1 = 10.25;
%   p.value2 = 2;
% 
%   p_int = double2int(p);
%   
%   p_int.value1
%       ans = 
%         10.2500
%
%   p_int.value2
%       ans = 
%         uint32
%         2
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

function [ p ] = double2int( p )
import utils.struc2cell

% convert structure to cell to get all fieldnames
fn = struc2cell(p);


for ii=1:length(fn)
    grps = strsplit(fn{ii}, '.');
    data_val = getfield(p, grps{:});
    
    if isnumeric(data_val) && isreal(data_val)
        try
            % check if data can converted to unsigned int
            if all(mod(data_val(:),1)==0) && all(all(data_val>=0))
                p = setfield(p,grps{:}, uint32(data_val));
             % if data is negative, convert it to int
            elseif all(mod(data_val(:),1)==0)
                p = setfield(p,grps{:}, int32(data_val));
            end
        catch
            keyboard
        end
    end
        
        
end
    
    
end







