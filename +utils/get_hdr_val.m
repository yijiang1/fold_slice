% Call function without arguments for instructions on how to use it

% Filename: $RCSfile: get_hdr_val.m,v $
%
% $Revision: 1.3 $  $Date: 2008/08/28 18:47:31 $
% $Author:  $
% $Tag: $
%
% Description:
% Find text signature in a bunch of cell strings from a file header and
% return the following value in the specified format. Example: 
% no_of_bin_bytes = get_hdr_val(header,'X-Binary-Size:','%f',1);
% The last parameter specifies whether the macro should exit with an error
% message if the text signature has not been found. 
%
% Note:
% Call without arguments for a brief help text.
%
% Dependencies:
% none
%
%
% history:
%
% May 7th 2008: 
% add number of input argument check and brief help text
%
% April 25th 2008: 1st version

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

function [outval,line_number,err] = get_hdr_val(header,signature,format,...
    exit_if_not_found)

% initialize output arguments
outval = 0;
line_number = 0;
err = 0;

if (nargin ~= 4)
    fprintf('Usage:\n');
    fprintf('[value,line_number,error]=%s(header,signature,format,exit_if_not_found);\n',...
        mfilename);
    fprintf('header             cell array with text lines as returned by cbfread or ebfread\n');
    fprintf('signature          string to be searched for in the header\n');
    fprintf('format             printf-like format specifier for the interpretation of the value that follows the signature\n');
    fprintf('exit_if_not_found  exit with an error in case either the signature or the value have not been found\n');
    error('Wrong number of input arguments.\n');
end

% search for the signature string
pos_found = strfind(header,signature);

% for sscanf the percentage sign has a special meaning
signature_sscanf = strrep(signature,'%','%%');

% loop over the search results for all header lines
for (ind=1:length(pos_found))
    % if the signature string has been found in this line
    if (length(pos_found{ind}) > 0)
        % get the following value in the specified format
        [outval,count] = sscanf(header{ind}(pos_found{ind}:end),...
            [signature_sscanf format]);
        % return an error if the signature and value combination has not
        % been found (i.e., the format specification did not match)
        if (count < 1)
            outval = 0;
            err = 1;
        else
            % return the first occurrence if more than one has been found
            if (count > 1)
                outval = outval(1);
            end
            % return the line number
            line_number = ind;
            return;
        end
    end
end

% no occurrence found
err = 1;
if (exit_if_not_found)
    error(['no header line with signature ''' signature ''' and format ' ...
        format ' found']);
end

return;
