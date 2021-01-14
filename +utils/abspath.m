% ABSPATH translate special symbols such as ~, ../, ./, in path to the absolute
% path 
%
% filename_with_path = abspath(filename_with_path)
% Inputs: 
%    **filename_with_path   original path 
%  Outputs: 
%    **filename_with_path   corrected path without special symbols 

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


function filename_with_path = abspath(filename_with_path)
    if isunix
        % replace home path 
        if startsWith(filename_with_path, '~/')
            filename_with_path = replace(filename_with_path, '~/', [char(java.lang.System.getProperty('user.home')), '/']); 
        end
        % replace root 
        nsteps = numel(strfind(filename_with_path, '../'));
        if nsteps > 0 && startsWith(filename_with_path, '../')
            new_path = pwd; 
            for ii = 1:nsteps
                new_path = fileparts(new_path); % remove the last folder from the path
            end
            filename_with_path = replace(filename_with_path, repmat('../', [1 nsteps]), [new_path, '/']);
        end
        % replace current folder 
        if  startsWith(filename_with_path, './')
            filename_with_path = replace(filename_with_path, './', [pwd, '/']); 
        end
    elseif ispc
        filename_with_path = replace(filename_with_path, '/', '\');
        filename_with_path = replace(filename_with_path, '.\', [pwd, '\']); 
    end
end

