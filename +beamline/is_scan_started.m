% out = is_scan_started(specDatFile,scanno)
% Detect 'X# ' in spec file to detect the end of a scan
% From spec compile  post_scan.mac

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

function out = is_scan_started(specDatFile,scanno)

    specDatFile = beamline.find_specDatFile(specDatFile);
    
    cmd = sprintf('grep -n ''#S %d'' %s', scanno,specDatFile);
    [~,sysout] = system(cmd);
    arrout = regexp(sysout,'[:\n ]','split');
    indS = find(strcmp(regexp(sysout,'[:\n ]','split'),'#S')==1); % Indices where #S is found
    for ii=indS % Loop over all #S found, this is to make sure we dont recognize S# 191 when looking for S# 19
        if strcmp(arrout(ii+1),sprintf('%d',scanno))
            out = true;
            return
        end
    end
    % If not found
    out = false;
    return
end
