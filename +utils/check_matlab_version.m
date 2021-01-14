%CHECK_MATLAB_VERSION check matlab version to make sure it is compatible
% ver...    compatible version number
%
% EXAMPLE:
%    check_matlab_version(9.2)
%
% MATLAB 9.0 - 2016a
% MATLAB 9.1 - 2016b
% MATLAB 9.2 - 2017a
% MATLAB 9.3 - 2017b


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

function check_matlab_version( ver )

current_version = version;
ver_str = strsplit(current_version, '.');
ver_num = str2double([ver_str{1} '.' ver_str{2} ver_str{3}]);

if ver_num < ver
    warning('You are using Maltab version %0.2f but the code was designed and tested with %0.2f.', ver_num, ver);
end

end

