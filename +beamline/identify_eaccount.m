% Return the current e-account user name in case this function is executed
% at the X12SA beamline, [] otherwise.

% Filename: $RCSfile: identify_eaccount.m,v $
%
% $Revision: 1.1 $  $Date: 2010/04/28 18:00:56 $
% $Author:  $
% $Tag: $
%
% Description:
% Return the current e-account user name in case this function is executed
% at the X12SA beamline, [] otherwise. 
%
% Note:
% none
%
% Dependencies:
% identify_system.m
%
%
% history:
%
% April 28th, 2010: 1st version

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

function [return_user_name] = identify_eaccount()
import utils.identify_system

persistent user_name;

if (isempty(user_name))
    user_name = [];
    % at the cSAXS beamline return the name of the current user as
    % e-account name
    sys_id = identify_system();
    if (strcmp(sys_id,'X12SA'))
        [st,un] = system('echo $USER');
        if ((st == 0) && (length(un) > 1))
            user_name = un(1:end-1);
        end
    end
end

return_user_name = user_name;
