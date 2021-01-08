% Call function without arguments for a detailed explanation of its use

% Filename: $RCSfile: datread.m,v $
%
% $Revision: 1.2 $  $Date: 2009/02/20 19:33:01 $
% $Author:  $
% $Tag: $
%
% Description:
% Macro for reading .dat files in self-defined data formats
%
% Note:
% Call without arguments for a brief help text.
%
% Dependencies:
% - image_read_set_default
% - fopen_until_exists
% - get_hdr_val
% - compiling cbf_uncompress.c increases speed but is not mandatory
%
%
% history:
%
% February 18th 2009: 1st version

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

function [frame,vararg_remain] = datread(filename,varargin)
import io.*
import utils.fopen_until_exists
import utils.get_hdr_val
import utils.char_to_cellstr

% 0: no debug information
% 1: some feedback
% 2: a lot of information
debug_level = 0;

% initialize return argument
frame = struct('header',[], 'data',[]);


% check minimum number of input arguments
if (nargin < 1)
    image_read_sub_help(mfilename,'cbf');
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
    error('The optional parameters have to be specified as ''name'',value pairs');
end
    
% set default values for the variable input arguments and parse the named
% parameters: 
vararg = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        otherwise
            % pass further arguments on to fopen_until_exists
            vararg{end+1} = name;
            vararg{end+1} = value;
    end
end


% expected maximum length for the text header
max_header_length = 4096;

% end-of-header signature
eoh_signature = [ '# end-of-header' char(10) ];

% try to open the data file
if (debug_level >= 1)
    fprintf('Opening %s.\n',filename);
end
[fid,vararg_remain] = fopen_until_exists(filename,vararg);
if (fid < 0)
    return;
end

% read all data at once
[fdat,fcount] = fread(fid,'uint8=>uint8');

% close input data file
fclose(fid);
if (debug_level >= 2)
    fprintf('%d data bytes read\n',fcount);
end

% files with header start with a has sign, otherwise just a series of
% numbers is expected
if (~strcmp(fdat(1),'#'))
    % convert the array to cell strings
    cell_values = char_to_cellstr( char(fdat)' );
    % convert the strings to double precision values
    frame.data = zeros(length(cell_values),1);
    for (ind = 1:length(cell_values))
        frame.data(ind) = str2double(cell_values{ind});
    end
    
    % No header information are available. 
    % Fake exposure time information to avoid problems in other
    % macros. 
    frame.header{end+1} = 'Exposure_time 1.0';
    % add the file modification date to the header
    dir_entry = dir(filename);
    frame.header{end+1} = [ 'DateTime ' dir_entry.date ];
end


% search for end of header signature within the expected maximum length of
% a header
end_of_header_pos = ...
    strfind( fdat(1:min(max_header_length,length(fdat)))',...
             eoh_signature );
if (length(end_of_header_pos) < 1)
    error( [ filename,': no header end signature found' ] );
    return;
end
if (debug_level >= 2)
    fprintf('Header length is %d bytes.\n',end_of_header_pos -1);
end

% return the complete header as lines of a cell array
frame.header = char_to_cellstr( char(fdat(1:(end_of_header_pos-1))') );

% increase the index to the first data byte
end_of_header_pos = end_of_header_pos + length(eoh_signature);

% check for information on the various dimensions in ascending speed order
dim1 = get_hdr_val(frame.header,'dim2','%d',1);
dim2 = get_hdr_val(frame.header,'dim1','%d',1);
dim3 = get_hdr_val(frame.header,'number-of-exposures','%d',1);
dim4 = get_hdr_val(frame.header,'channels','%d',1);

if (debug_level >= 2)
    fprintf('Frame dimensions are %d x %d % %d x %d.\n', ...
            dim4,dim3,dim2,dim1);
end

% store the numbers in the array, fastest axis first
frame.data = zeros(dim4,dim3,dim2,dim1);

% convert the strings to double precision values
data_1d = sscanf(char(fdat(end_of_header_pos:end))','%f');
frame.no_of_el_read = length(data_1d);
if (length(data_1d) > numel(frame.data))
    frame.data = zeros(dim4,dim3,dim2,ceil(length(data_1d)/(dim4*dim3*dim2)));
end
frame.data(1:frame.no_of_el_read) = data_1d;

