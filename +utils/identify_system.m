% Identify the current system to set useful default parameter values in
% default_parameter_value.m. 
% A modified version of both macros at the beginning of the Matlab search
% path may be used to define local standard parameters. 

% Filename: $RCSfile: identify_system.m,v $
%
% $Revision: 1.3 $  $Date: 2010/07/22 15:08:21 $
% $Author:  $
% $Tag: $
%
% Description:
% Identify the current system to set useful default parameter values in
% default_parameter_value.m. 
% A modified version of both macros at the beginning of the Matlab search
% path may be used to define local standard parameters. 
%
% Note:
% none
%
% Dependencies:
% none
%
%
% history:
%
% June 2nd 2009:
% buffer current system ID for later calls to speed up execution
%
% April 16th, 2009: 1st version

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

function [return_system_id_str return_other_system_flags] = identify_system()

persistent system_id_str;
persistent parallel_computing_toolbox_available;

if (isempty(system_id_str))
    % default value
    system_id_str = 'other';

    if (isunix)
        % check for a known network name of the PC Matlab is running on
        [status,hostname] = unix('hostname');
        if (status == 0)
            hostname = sscanf(hostname,'%s');
            if length(hostname)>4 && strcmp(hostname(1:5),'x12sa')
                system_id_str = 'X12SA';
            else
                switch hostname 
                    case {'pc6024', 'pc5369'}
                        system_id_str = 'DPC lab';
                    case {'mpc1054'}
                        system_id_str = 'mDPC lab';
                    case {'pc5211', 'mpc1144', 'mpc1145'}
                        system_id_str = 'cSAXS-mobile';
                    case {'lccxs01', 'lccxs02', 'mpc1208'}
                        system_id_str = 'CXS compute node';
                end
            end
        end
    else
        % neither Linux nor Mac
        system_id_str = 'Windows';
    end
end

% check for the parallel computing toolbox being available
if (isempty(parallel_computing_toolbox_available))
    parallel_computing_toolbox_available = false;
    
    versions = ver;
    for line = 1:length(versions)
        if strfind(versions(line).Name, 'Parallel Computing Toolbox')
            parallel_computing_toolbox_available = true;
        end
    end
end


% compile return values

return_system_id_str = system_id_str;

return_other_system_flags.parallel_computing_toolbox_available = parallel_computing_toolbox_available;
