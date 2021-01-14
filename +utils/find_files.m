% Call function without arguments for a detailed explanation of its use

% Filename: $RCSfile: find_files.m,v $
%
% $Revision: 1.9 $  $Date: 2012/08/07 16:39:30 $
% $Author:  $
% $Tag: $
%
% Description:
% find file names matching the specified mask
%
% Note:
% Call without arguments for a brief help text.
%
% Dependencies: 
% - Linux/Unix find command, if specified to use
%
% history:
%
% October 10th 2009:
% only check for files if changing to the directory was possible
%
% September 14th 2008: 
% bug fix: add directory to filename in isdir check
%
% September 4th 2008:
% bug fix: vararg_remain was not filled and unhandled parameters did not
% cause an error
%
% June 16th 2008: send find output through sort
%
% June 10th 2008: 1st documented version

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

function [ directory, fnames, vararg_remain ] = find_files( filename_mask, varargin )
import io.image_read
import utils.default_parameter_value

% initialize return arguments
fnames = [ ];

% set default values
use_find = default_parameter_value(mfilename,'UseFind');
unhandled_par_error = default_parameter_value(mfilename,'UnhandledParError');

% check minimum number of input arguments
if (nargin < 1)
    fprintf('\nUsage:\n');
    fprintf('[directory filenames]=%s(filename_mask,  [[,<name>,<value>] ...]);\n',mfilename);
    fprintf('filename_mask can be something like ''*.cbf'' or ''image.cbf''\n');
    fprintf('The optional <name>,<value> pairs are:\n');
    fprintf('''UseFind'',<0-no, 1-yes>              use Linux/Unix command find to interprete the filename mask, default is %d\n',use_find);
    fprintf('''UnhandledParError'',<0-no,1-yes>     exit in case not all named parameters are used/known, default is %d\n',unhandled_par_error);
    fprintf('Examples:\n');
    fprintf('%s(''~/Data10/pilatus/mydatadir/*.cbf'',''OutdirData'',''~/Data10/analysis/my_int_dir/'');\n',mfilename);
    fprintf('Additional <name>,<value> pairs recognized by image_read can be specified.\n');
    error('At least the filename mask has to be specified as input argument.');
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

% parse the variable input arguments:
% initialize the list of unhandled parameters
vararg_remain = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'UseFind' 
            use_find = value;
        case 'UnhandledParError'
            unhandled_par_error = value;
        otherwise
            vararg_remain{end+1} = name; %#ok<AGROW>
            vararg_remain{end+1} = value; %#ok<AGROW>
    end
end

[directory, name, ext] = fileparts(filename_mask);
% add slash to directories
if ((~isempty(directory)) && (directory(end) ~= '/'))
    directory = [ directory '/' ];
end

% exit in case of unhandled named parameters, if this has not been switched
% off
if ((unhandled_par_error) && (~isempty(vararg_remain)))
    vararg_remain %#ok<NOPRT>
    error('Not all named parameters have been handled.');
end


% search matching filenames
if (use_find)
    find_cmd = sprintf('find . -noleaf -maxdepth 1 -name ''%s''|sort',[ name ext ]);
    cd_cmd = '';
    if (~isempty(directory))
        %Note by YJ: different linux accounts use different "cd" commands.
        %cd_cmd may cause error for some users (e.g. user2idd)
        cd_cmd = sprintf('cd %s 2>/dev/null',directory);
    end
    
    % if the directory exists check for files within it
    st = 1;
    if ((isempty(directory)) || (exist(directory,'dir')))
        [st,files]=system([cd_cmd ';' find_cmd ]);
        %disp([cd_cmd ';' find_cmd ])
        %disp(files)
    end
    % store names of files found in fnames
    if (st == 0)
        % count number of newline characters
        no_of_files = length(sscanf(files,'%*[^\n]%1c'));
        fnames = struct('name',cell(1,no_of_files),'isdir',cell(1,no_of_files));
        % extract file names
        file_ind = 1;
        while (~isempty(files))
            name = sscanf(files,'%[^\n]',1);
            files = files( (length(name)+2):end );
            if ((length(name) > 2) && (strcmp(name(1:2),'./')))
                name = name(3:end);
            end
            % exclude directory entries . and .. and error messages
            % starting with find: that may occur if temporary Pilatus files
            % vanish
            if ((~strcmp(name,'.')) && ...
                ((length(name) < 5) || (~strcmp(name(1:5),'find:'))))
                fnames(file_ind).name = name;
                fnames(file_ind).isdir = isfolder([ directory fnames(file_ind).name ]);
                file_ind = file_ind +1;
            end
        end
        % shorten the result if for example '.' entries have been skipped
        if (file_ind <= no_of_files)
            fnames = fnames(1:(file_ind-1));
        end
    end
else
    fnames = dir( filename_mask );
end
