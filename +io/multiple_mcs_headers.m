% [] = multiple_mcs_headers(scan_no_from,scan_no_to)

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

function [] = multiple_mcs_headers(scan_no_from,scan_no_to)
import utils.compile_x12sa_filename

if (nargin ~= 2)
    fprintf('Usage: %s <scan no. from> <scan no. to>\n',mfilename);
    fprintf('Searches for multiple MCS headers, renames the matching files and\n');
    fprintf('writes a new file with only the last header and its data.\n');
    fprintf('This is just a workaround for the current implementation of the MCS\n');
    fprintf('where spec stores the MCS data several times in case a cont_line is repeated.\n');
end

for (scan_no = scan_no_from:scan_no_to)
    dirname = compile_x12sa_filename(scan_no, -1, 'BasePath','~/Data10/mcs/');
    dirinfo = dir([ dirname '*.dat' ]);
    if (length(dirinfo) < 1)
        fprintf('%s not found\n',filename);
    else
        for (ind = 1:length(dirinfo))
            filename = [dirname dirinfo(ind).name];
            fid = fopen(filename,'r');
            % read all data at once
            [fdat,~] = fread(fid,'uint8=>uint8');
            fclose(fid);
            % check for multiple fileheaders
            header_start = strfind(fdat','# MCS file version');
            if (length(header_start) > 1)
                fprintf('%s has %d headers\n',filename, length(header_start));
                % rename the file
                movefile(filename, [filename '_org']);
                % use the last file part
                fdat = fdat(header_start(end):end);
                % store the truncated data
                fid = fopen(filename,'w');
                fwrite(fid,fdat);
                fclose(fid);
            end
        end
    end
end
