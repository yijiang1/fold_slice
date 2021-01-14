% Call function without arguments for a detailed explanation of its use

% Filename: $RCSfile: fopen_until_exists.m,v $
%
% $Revision: 1.9 $  $Date: 2011/08/13 17:37:15 $
% $Author:  $
% $Tag: $
%
% Description:
% Open a file, in case of failure retry repeatedly if this has been
% specified. 
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
% September 5th 2009:
% bug fix in the zero file length check
%
% August 28th 2008: 
% use dir rather than fopen to check for the file and check additionally
% that it is not of length zero
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

function [fid,vararg_remain] = fopen_until_exists(filename,varargin)

% set default values for the variable input arguments and parse the named
% parameters:

% If the file has not been found and if this value is greater than 0.0 than
% sleep for the specified time in seconds and retry reading the file. 
% This is repeated until the file has been successfully read
% (retry_read_max=0) or until the maximum number of iterations is exceeded
% (retry_read_max>0). 
retry_read_sleep_sec = 0.0;
retry_read_max = 0;
retry_sleep_when_found_sec = 0.0;

% exit with error message if the file has not been found
error_if_not_found = 1;

% display a message once in case opening failed
message_if_not_found = 1;

if (nargin < 1)
    fprintf('Usage:\n');
    fprintf('[fid] = %s(filename [[,<name>,<value>],...]);\n',...
        mfilename);
    fprintf('filename                             name of the file to open\n');
    fprintf('The optional name value pairs are:\n');
    fprintf('''RetryReadSleep'',<seconds>           if greater than zero retry opening after this time (default: 0.0)\n');
    fprintf('''RetryReadMax'',<0-...>               maximum no. of retries, 0 for infinity (default: 0)\n');
    fprintf('''RetrySleepWhenFound'',<seconds>      if greater than zero wait for this time after a retry succeeded (default: %.1f)\n', ...
        retry_sleep_when_found_sec);
    fprintf('''MessageIfNotFound'',<0-no,1-yes>     display a mesage if not found, 1-yes is default\n');
    fprintf('''ErrorIfNotFound'',<0-no,1-yes>       exit with an error if not found, default is 1-yes\n');
    fprintf('The file ID of the opened file is returned or -1 in case of failure.\n');
    error('Invalid number of input parameters.');
end

% check minimum number of input arguments
if (nargin < 1)
    display_help();
    error('At least the filename has to be specified as input parameter.');
end

% accept cell array with name/value pairs as well
no_of_in_arg = nargin;
if (nargin == 2)
    if (isempty(varargin))
        % ignore empty cell array
        no_of_in_arg = no_of_in_arg -1;
    else
        if (iscell(varargin{1}))
            % use a filled one given as first and only variable parameter
            varargin = varargin{1};
            no_of_in_arg = 1 + length(varargin);
        end
    end
end
% check number of input arguments
if (rem(no_of_in_arg,2) ~= 1)
    error('The optional parameters have to be specified as ''name'',''value'' pairs');
end


% parse the variable input arguments
vararg_remain = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'RetryReadSleep'
            retry_read_sleep_sec = value;
        case 'RetryReadMax'
            retry_read_max = value;
        case 'RetrySleepWhenFound'
            retry_sleep_when_found_sec = value;
        case 'MessageIfNotFound'
            message_if_not_found = value;
        case 'ErrorIfNotFound'
            error_if_not_found = value;
        otherwise
            vararg_remain{end+1} = name;
            vararg_remain{end+1} = value;            
    end
end


% try to access the file entry
file_non_empty = 0;
dir_entry = dir(filename);

% if it has not been found or if it is empty
if ((isempty(dir_entry)) || (size(dir_entry,1) == 0) || ...
    (dir_entry.bytes <= 0))
    if (message_if_not_found)
        if (isempty(dir_entry))
            fprintf('%s not found',filename);
        else
            fprintf('%s found but of zero length',filename);            
        end
    end
    % retry, if this has been specified
    if (retry_read_sleep_sec > 0.0)
        if (message_if_not_found)
            fprintf(', retrying\n');
        end
        % repeat until found or the specified number of repeats has been
        % exceeded (zero repeats means repeat endlessly)
        retry_read_ct = retry_read_max;
        while ((~file_non_empty) && ...
               ((retry_read_max <= 0) || (retry_read_ct > 0)))
            fprintf('Pausing %d seconds and retrying \n',retry_read_sleep_sec);
            pause(retry_read_sleep_sec);
            dir_entry = dir(filename);
            if ((~isempty(dir_entry)) && (dir_entry.bytes > 0))
                file_non_empty = 1;
                % workaround option for various problems, 
                % not for permanent use
                if (retry_sleep_when_found_sec > 0)
                    pause(retry_sleep_when_found_sec);
                end
            end
            retry_read_ct = retry_read_ct -1;
        end
    else
        fprintf('\n');
    end
else
    file_non_empty = 1;
end

% open the file for read access
if (file_non_empty)
    fid = fopen(filename,'r');
else
  fid = -1;
end

% exit with an error message, if this has been specified and if the file
% could not be opened 
if (fid < 0)
    if (error_if_not_found)
        ME = MException('fopen_until_exists:not_found', ...
            strjoin({'File', filename, 'does not exist'}));
        throwAsCaller(ME);
    end
end
