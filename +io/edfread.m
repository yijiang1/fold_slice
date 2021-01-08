% Call function without arguments for a detailed explanation of its use

% Filename: $RCSfile: edfread.m,v $
%
% $Revision: 1.1 $  $Date: 2008/06/10 17:05:14 $
% $Author:  $
% $Tag: $
%
% Description:
% Macro for reading ESRF Data Format (EDF) files written by the
% Pilatus detector control program camserver. 
%
% Note:
% Currently this routine supports only the subset of EDF features needed to
% read the Pilatus detector data. 
% Call without arguments for a brief help text.
%
% Dependencies:
% - fopen_until_exists
% - get_hdr_val
% - image_read_set_default
%
%
% history:
%
% May 9th 2008: 1st version after redesign

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

function [frame,vararg_remain] = edfread(filename,varargin)
import io.*
import utils.char_to_cellstr
import utils.fopen_until_exists
import utils.get_hdr_val

% 0: no debug information
% 1: some feedback
% 2: a lot of information
debug_level = 0;

% initialize return argument
frame = struct('header',[], 'data',[]);

% check minimum number of input arguments
if (nargin < 1)
    image_read_sub_help(mfilename,'edf');
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

% search for end of header signature within the expected maximum length of
% a header
end_of_header_pos = 1024;
max_pos = min( max_header_length,length(fdat) );

while ((end_of_header_pos < max_pos) && (fdat(end_of_header_pos-1) ~= '}'))
    end_of_header_pos = end_of_header_pos +1024;
end
if (end_of_header_pos >= max_pos)
    error('no header end signature found');
end
if (debug_level >= 2)
    fprintf('Header length is %d bytes.\n',end_of_header_pos);
end
data_length = fcount - end_of_header_pos;

% convert the header to lines of a cell array
frame.header = char_to_cellstr( char(fdat(1:(end_of_header_pos-1))') );

% check for opening parenthesis
if (frame.header{1} ~= '{')
    error([filename ': EDF start ,''{'' not found in first line ''' ...
        frame.header{1} '''' ]);
end


% extract the mandatory information for data extraction from the header:
byte_order = get_hdr_val(frame.header,'ByteOrder',' = %s',1);
dim1 = get_hdr_val(frame.header,'Dim_1',' = %d',1);
dim2 = get_hdr_val(frame.header,'Dim_2',' = %d',1);
data_type = get_hdr_val(frame.header,'DataType',' = %s',1);
if (debug_level >= 2)
    fprintf('Byte order is %s\n',byte_order);
    fprintf('Frame dimensions are %d x %d.\n',dim2,dim1);
    fprintf('Data type is %s\n',data_type);
end

% determine number of bytes per pixel
switch data_type
    case 'UnsignedByte',
        bytes_per_pixel = 1;
        data_class = 'uint8';
    case 'UnsignedShort',
        bytes_per_pixel = 2;
        data_class = 'uint16';
    case {'SignedInteger','UnsignedInteger','UnsignedInt','UnsignedLong'}
        bytes_per_pixel = 4;
        data_class = 'uint32';
    case {'Float','FloatValue','Real'}
        bytes_per_pixel = 4;
        data_class = 'single';
    case 'DoubleValue'
        bytes_per_pixel = 8;
        data_class = 'double';
    otherwise
        error('unsupported data type %s',data_type);
end
no_of_bytes = bytes_per_pixel * dim1 * dim2;
if (debug_level >= 2)
    fprintf('%d bytes per pixel, %d in total expected, %d available\n',...
        bytes_per_pixel,no_of_bytes,data_length);
end

% check length of available data
if (no_of_bytes > data_length)
    error('%d data bytes expected, %d are available',...
        no_of_bytes,data_length);
end

% compare file with machine byte order, swap if necessary
[str,maxsize,endian] = computer;
if (((strcmp(byte_order,'HighByteFirst')) && (endian == 'L')) || ...
    ((strcmp(byte_order,'LowByteFirst')) && (endian == 'H')))
    if (debug_level >= 2)
        fprintf('Machine byte order is %s: swapping data bytes\n',...
            endian,bytes_per_pixel);
    end
    dat = fdat(end_of_header_pos+1:end_of_header_pos+no_of_bytes);
    dat = reshape(dat,bytes_per_pixel,[]);
    dat = flipud(dat);
    fdat(end_of_header_pos+1:end_of_header_pos+no_of_bytes) = dat(:);
end

% extract the frame from the binary data
[frame.data] = ...
    double(reshape(typecast(fdat(end_of_header_pos+1:end_of_header_pos+no_of_bytes),...
                            data_class),...
                   dim1,dim2));

    
% conversion to standard view on Pilatus 2M data at the SLS/cSAXS beamline
% if (~original_orientation)
%     % this is slow, even slower is fliplr(flipud(frame.'))
%     frame.data = frame.data(end:-1:1,end:-1:1)';
% end
