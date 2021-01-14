
% FUNCTION [mem_avail, mem_total] = check_availible_memory()
% get availible free memory in linux in MB 
    
    
%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2018 by Paul Scherrer Institute (http://www.psi.ch)    |
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


function varargout =check_available_memory()

    if isunix
        meminfo = importdata('/proc/meminfo'); 
        mem_total = str2num((regexprep(meminfo.textdata{1}, '[a-zA-Z \:]', '')))/1e3; 
        mem_avail = str2num((regexprep(meminfo.textdata{3}, '[a-zA-Z \:]', '')))/1e3; 


        if nargout > 0
            vlevel = 2;
        else
            vlevel = 0; 
        end
        utils.verbose(vlevel, '====   %.0f GB == %.0f%% RAM free ====', mem_avail/1e3, mem_avail/mem_total*100);

        if mem_avail/mem_total < 0.2
            warning('Less than 20% RAM left..');
            !free -h    
        end
    elseif ispc
        [~,sV] = memory;
        mem_avail = sV.PhysicalMemory.Available/1e6; 
        mem_total = sV.PhysicalMemory.Total/1e6; 
    else
        error('Unsupported architecture')
    end
    
    
    if nargout > 0
        varargout = {mem_avail, mem_total};
    end


end
