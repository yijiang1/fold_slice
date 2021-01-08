% [] = image_read_sub_help(m_file_name,extension,varargin)
% parameter help for sub-routines of image_read like cbfread, edfread,
% speread and fliread

% Filename: $RCSfile: image_read_sub_help.m,v $
%
% $Revision: 1.1 $  $Date: 2008/06/10 17:05:14 $
% $Author:  $
% $Tag: $
%
% Description:
% parameter help for sub-routines of image_read like cbfread, edfread,
% speread and fliread
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
% May 9th 2008: 1st version

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

function [] = image_read_sub_help(m_file_name,extension,varargin)
import io.*
% check minimum number of input arguments
if (nargin < 2)
    error('At least the m-file name and the extension have to be specified as input parameter.');
end

% accept cell array with name/value pairs as well
no_of_in_arg = nargin;
if (nargin == 3)
    if (isempty(varargin))
        % ignore empty cell array
        no_of_in_arg = no_of_in_arg -1;
    else
        if (iscell(varargin{1}))
            % use a filled one given as first and only variable parameter
            varargin = varargin{1};
            no_of_in_arg = 2 + length(varargin);
        end
    end
end

% check number of input arguments
if (rem(no_of_in_arg,2) ~= 0)
    error('The optional parameters have to be specified as ''name'',''value'' pairs');
end


% parse the variable input arguments
examples = 1;
extension_returned = 0;
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'Examples' 
            examples = value;
        case 'ExtensionReturned'
            extension_returned = 1;
        otherwise
            error('unknown parameter %s',name);
    end
end

fprintf('Usage:\n')
if (extension_returned)
    fprintf('[frame]=%s(<filename> [[,<name>,<value>] ...]);\n',...
        m_file_name);
else
    fprintf('[frame]=%s(<filename> [[,<name>,<value>] ...]);\n',...
        m_file_name);
end
fprintf('The optional <name>,<value> pairs are:\n');
fprintf('''RetryReadSleep'',<seconds>           if greater than zero retry opening after this time (default: 0.0)\n');
fprintf('''RetryReadMax'',<0-...>               maximum no. of retries, 0 for infinity (default: 0)\n');
fprintf('''MessageIfNotFound'',<0-no,1-yes>     display a mesage if not found, 1-yes is default\n');
fprintf('''ErrorIfNotFound'',<0-no,1-yes>       exit with an error if not found, default is 1-yes\n');
if (examples)
    fprintf('\n');
    fprintf('Examples:\n');
    fprintf('[frame]=%s(''~/Data10/roper/image.%s'');\n',...
        m_file_name,extension);
    fprintf('\n');
    fprintf('The returned structure has the fields data, header and extension.\n');
end
