% REPORT_GPU_USAGE check who is using selected GPU and how much memory the
% user uses 
%
% report_GPU_usage(gpu_id)
% 
% Inputs: 
%   **gpu_id - ID of the checked GPU

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



function report_GPU_usage(gpu_id)
    
    if nargin == 0
        gpu = gpuDevice; 
    else
        gpu = gpuDevice(gpu_id);
    end

    % find user to whom belong processes on the selected GPU 
    [~,out] = system(sprintf('nvidia-smi | awk ''$2=="Processes:" {p=1} p && $2 == %i && $3 > 0 {print $3}''',gpu.Index-1)); 
    pids = str2num(out);
    % find memory use by the processes %     polyfit_order - -1 = dont assume anything about the removed phase, remove line fit in every horizontal line separatelly  
    [~,out] = system(sprintf('nvidia-smi | awk ''$2=="Processes:" {p=1} p && $2 == %i && $3 > 0 {print $6}''',gpu.Index-1));
    out = splitlines(out); 
    for ii = 1:length(out)-1
        memory_use(ii) = str2num(replace(out{ii}, 'MiB', ''));
    end
    
    for ii = 1:length(pids)
      [~,out] = system(sprintf('ps aux | grep %i | awk ''{print $1}''| head -n1', pids(ii)));
       user_list{ii} = out(1:end-1);
    end
    users = unique(user_list); 
    for ii = 1:length(users)
        user_memory(ii) = sum(memory_use(ismember(user_list, users{ii}))); 
    end
    
    utils.verbose(0,'Following users have processes on the selected GPU %i:', gpu.Index);
    for  ii = 1:length(users)
        utils.verbose(0, 'User: %s \t Used memory: %4.3gGB', users{ii} ,user_memory(ii)/1e3)
    end
    
end