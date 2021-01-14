% CHECK_CPU_LOAD returns user cpu usage of specified hosts
% 
% hosts (optional)...       list of nodes; use 'x12sa' for all x12sa nodes
% used_nodes (optional)...  prints warning/summary for used nodes (default: true)
% ssh_auth (optional) ...   system echo if ssh authentication is needed (default: false)
% thr (optional)...         set threshold for used nodes (default: 15 %)
% vm_cycles (optional)...   number of cycles for cpu usage (default: 2)

% 03/2017

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

function [ cpu_load_bl, any_used_cpu ] = check_cpu_load( varargin)

any_used_cpu = false;

[~, hostname] = system('hostname');
host_pre = strsplit(hostname, '-');
switch host_pre{1}
    case 'ra'
        vmstat_nr = 13;
    case 'x12sa'
        vmstat_nr = 15;
    otherwise
        vmstat_nr = 15;
end

if nargin < 1

    switch host_pre{1}
        case 'x12sa'
            hosts = {'x12sa-cn-1', 'x12sa-cn-2', 'x12sa-cn-3', 'x12sa-cn-4', 'x12sa-cn-5', 'x12sa-cn-6'};
        otherwise
            error('Please specify your hosts.');
    end
else
    if ischar(varargin{1})
        varargin{1}  = strcell(varargin{1});  % backward compatibility 
    end
    hosts = varargin{1};
    if strcmp(varargin{1}, 'x12sa')
        hosts = {'x12sa-cn-1', 'x12sa-cn-2', 'x12sa-cn-3', 'x12sa-cn-4', 'x12sa-cn-5', 'x12sa-cn-6'};
    end
end


% check if warnings/summary are needed
if nargin > 1
    used_nodes = varargin{2};
else
    used_nodes = true;
end

% check ssh_auth
if nargin > 2
    ssh_auth = varargin{3};
else
    ssh_auth = false;
end

% check if thr for cpu load is specified
if nargin > 3
    thr = varargin{4};
else
    thr = 15;
end

% check if cycles are specified
if nargin > 4
    top_cycles = varargin{5};
else
    top_cycles = 2;
end

cycles = sprintf('%i', top_cycles);

cpu_load_bl = zeros(length(hosts),1);

% ssh to hosts and check cpu load
for i=1:size(hosts,1)
    
    ssh_call_cpu = ['ssh ' hosts{i}, ' vmstat 1 ' cycles ' | tail -1 | awk ''{print 100 - $' num2str(vmstat_nr) '}'''];

    if ssh_auth
        [~, result] = system(ssh_call_cpu, '-echo');
    else
        [~, result] = system(ssh_call_cpu);
    end
    res = strsplit(result, '\n');
    for j=1:size(res,2)
        if ~isnan(str2double(res{j}))
            cpu_load_bl(i) = str2double(res{j});
            if used_nodes
                fprintf('%s:  %i%%\n', hosts{i}, cpu_load_bl(i));
            end
        end
    end
end


% print warnings and summary if needed
if used_nodes
    for i=1:size(hosts)
        if cpu_load_bl(i) > thr
            fprintf('Host %s is currently used!\n', hosts{i});
            any_used_cpu = true;
        end
    end
end

end

