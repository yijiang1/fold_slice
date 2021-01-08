% GET_ROI find optimal rectange containing the region given by a binary mask 
%
%  [ROI] = get_ROI(mask, extent, multiply_of)
%
% Inputs: 
%   **mask      - binary mask, true for valid region 
%   **extent    - extra region around the mask calculated as extent*[W,H], default == 0
%   **multiply_of - make the selected ROI dividable by this number, default = 1
% *returns*: 
%   ++ROI       - 2x1 cell containing indices of the selected region 


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

function [ROI] = get_ROI(mask, extent, multiply_of)


    if all(mask(:) == 0)
        error('Empty mask')
    end

    if nargin < 3
        multiply_of = 1;
    end
    if nargin < 2
       extent = 0;
    end
    
    if ~islogical(mask)
       error('Not implemented')
    end

    x = any(mask,2);
    y = any(mask,1);
    coord = gather([find(x, 1,'first'), find(x, 1,'last'), find(y, 1,'first'), find(y, 1,'last')]);

    w = (coord(2) - coord(1));
    h = (coord(4) - coord(3));
    Cx = (coord(2) +coord(1))/2;
    Cy = (coord(4) + coord(3))/2;

    %YJ: HAS BUG - WRONG INDEX FOR ENTIRE FOV. 
    %CHEAP FIX: USE A SMALL NEGATIVE EXTENT
    
    %why only use extent for horitonal direction? 
    coord(1) =  floor(Cx  - ceil( (0.5 + extent) *w ));
    coord(2) =  ceil(Cx  + ceil((0.5 + extent) * w ));
    %% added by YJ to apply extend to vertical direction
    coord(3) =  floor(Cy  - ceil( (0.5 + extent) *h ));
    coord(4) =  ceil(Cy  + ceil((0.5 + extent) * h )) ;
    %%
    w = floor((coord(2) - coord(1))/multiply_of)*multiply_of;
    h = floor((coord(4) - coord(3))/multiply_of)*multiply_of;
    
    coord(2) = coord(1)+w-1;
    coord(4) = coord(3)+h-1;

    coord([1,3]) = max(1,coord([1,3]));
    coord([2,4]) = min(size(mask),coord([2,4]));

    ROI = {coord(1):coord(2), coord(3):coord(4)};

end