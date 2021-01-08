%FIND_BASE_PACKAGE
% finds the path to the cSAXS base package by looking for a specific file (+math)
% the code goes up to 3 levels down in the folder structure and it tries to find any folder matching ./*/+math/
%
% returns:
% ++ base_package_path              path to the cSAXS base package
%
% Example how to use it: addpath(find_base_package())


%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
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



function base_package_path = find_base_package()
    maxdepth = 3;
	test_path = '+math'; % one file to find them all
    
    lvl = 1;
    cpath = '';
    ret = '';
    while isempty(ret) && ~contains(strtrim(ret), test_path)
        [~, ret] = system(sprintf('find -L %s -maxdepth 2 -type d -name "%s"', cpath, test_path));
        if lvl > maxdepth
            break
        end
        lvl = lvl + 1;
        cpath = [cpath '../'];
    end
    ret = split(ret); 
    base_package_path = strtrim(ret{1}); 
    base_package_path = base_package_path(1:end-length(test_path));

    if isempty(base_package_path)
        error('cSAXS base package was not found')
    end
end
