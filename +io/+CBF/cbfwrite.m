% Call function without arguments for a detailed explanation of its use

% Filename: $RCSfile: cbfwrite.m,v $
%
% $Revision: 1.2 $  $Date: 2014/04/17 16:49:22 $
% $Author:  $
% $Tag: $
%
% Description:
% Macro for writing data to a Crystallographic Binary File (CBF) file. 
% Assumes that the data have been read from such a file beforehand, i.e.,
% that the header exists already. 
%
%
% history:
%
% April 10th 2014: 1st version

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

function [fcount_total] = cbfwrite(filename,frame,varargin)
import io.CBF.*
import io.*
import utils.default_parameter_value
import utils.fopen_until_exists
import utils.get_hdr_val

% 0: no debug information
% 1: some feedback
% 2: a lot of information
debug_level = 0;

% determine default orientation based on the file name extension
orient_by_extension = default_parameter_value('image_read','OrientByExtension');

fcount_total = -1;

% check minimum number of input arguments
if (nargin < 2)
    fprintf('\nUsage:\n');
    fprintf('[bytes_written]=%s( filename, frame-structure [[,<name>,<value>]...]);\n',mfilename);
    fprintf('Write a CBF file from data that have been read via cbfread.m before.\n');
    fprintf('Conventions are likely to be specific for the PILATUS detector.\n')
    fprintf('\n');
    fprintf('The optional <name>,<value> pairs are:\n');
    image_orient_help(mfilename,'ParametersOnly',1);
    fprintf('\n');
    fprintf('The file name should be the name of a single file without wildcards.\n');

    error('At least the filename and the data-structure have to be specified as input parameter.');
end

% check frame structure
if (~isstruct(frame))
    error('The 2nd parameter needs to be a structure as read by cbfread.m');
end
if (~isfield(frame,'header'))
    error('The 2nd parameter needs to be a structure with a field ''header'', as read by cbfread.m');
end
if (~isfield(frame,'data'))
    error('The 2nd parameter needs to be a structure with a field ''data'', as read by cbfread.m');
end
if (size(frame.header,2) ~= 1)
    error('The header has the wrong dimensions. Only single-frame structures are suppoerted.\n');
end

% check minimum number of input arguments
if (nargin < 2)
    error('At least the filename and the frame have to be specified as input parameter.');
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
    error('The optional parameters have to be specified as ''name'',value pairs');
end
    
% set default values for the variable input arguments and parse the named
% parameters: 
vararg_remain = cell(4,1);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'OrientByExtension' 
            orient_by_extension = value;
        otherwise
            % pass further arguments on to fopen_until_exists
            vararg_remain{end+1} = name; %#ok<AGROW>
            vararg_remain{end+1} = value; %#ok<AGROW>
    end
end
vararg_remain{1} = 'OrientByExtension';
vararg_remain{2} = orient_by_extension;
vararg_remain{3} = 'InvertOrientation';
vararg_remain{4} = 1;

% get the length of the zero-padding at the end of the file
padding_length = get_hdr_val(frame.header{1},'X-Binary-Size-Padding:','%d',1);

% end of header signature
eoh_signature = char([ 12 26 4 213 ]);

% orient image
[frame,vararg_remain] = image_orient(frame,vararg_remain);


% open file for write-access, overwrite in case of an existing file
[fid] = fopen(filename,'w');
if (fid < 0)
    return;
end


% In byte-offset compression the difference to the previous pixel value is
% stored as a byte, 16-bit integer or 32-bit integer, depending on its
% size. 
% The sizes above one byte are indicated by the escape sequence -1 in the
% previous data format, i.e, a 32-bit integer is preceded by the sequence %
% 0x80 (too large for a byte)
% 0x8000 (too large for a 16-bit integer). 
ind_out = 1;
no_of_pixels = size(frame.data,1) * size(frame.data,2);

% initialize output array at maximum size
frame_out = zeros(no_of_pixels*4 + padding_length,1,'uint8');
val_prev = 0;
for ind_in = 1:no_of_pixels
    val_diff = frame.data(ind_in) - val_prev;
    if (abs(val_diff) < 128)
        % write differences in the range from -127 to 127 directly as 8bit
        % signed integer:
        % manual complement to emulate the sign
        if (val_diff < 0)
            val_diff = 256 + val_diff;
        end
        frame_out(ind_out) = val_diff;
        ind_out = ind_out +1;
    else
        % escape with 0x80
        frame_out(ind_out) = 128;
        % check for 16-bit integer value
        if (abs(val_diff) < 32768)
            % write signed 16bit integer value:
            % manual complement to emulate the sign
            if (val_diff < 0)
                val_diff = 65536 + val_diff;
            end
            frame_out(ind_out+1) = bitand(val_diff,255);
            frame_out(ind_out+2) = bitand(val_diff,65280)/256;
            ind_out = ind_out +3;
        else
            % escape with 0x8000
            frame_out(ind_out+1) = 0;
            frame_out(ind_out+2) = 128;
            % write signed 32bit integer value:
            % manual complement to emulate the sign
            if (val_diff < 0)
                val_diff = 4294967296 + val_diff;
            end
            frame_out(ind_out+3) = bitand(val_diff,255);
            frame_out(ind_out+4) = bitand(val_diff,65280)/256;
            frame_out(ind_out+5) = bitand(val_diff,16711680)/65536;
            frame_out(ind_out+6) = bitand(val_diff,4278190080)/16777216;
            ind_out = ind_out +7;
        end
    end
    % the current intensity value becomes the previous value
    val_prev = frame.data(ind_in);
end

% calculate length of binary data including zero-padding at the end
binary_length = ind_out -1 + padding_length;

% update this value in the file-header
[~, line_no] = get_hdr_val(frame.header{1},'X-Binary-Size:','%f',1);
frame.header{1}{line_no} = sprintf('X-Binary-Size: %d',binary_length-padding_length);

% write header
[fcount_total] = fprintf(fid,'%s\r\n',frame.header{1}{:});
if (fcount_total == 0)
    error('Could not write any header data.\n');
end
% write end-of-header signature
[fcount] = fwrite(fid,eoh_signature,'uint8');
fcount_total = fcount_total + fcount;
if (fcount == 0)
    error('Could not write end-of-header signature.\n');
end

% write binary data
[fcount] = fwrite(fid,frame_out(1:(binary_length)));
fcount_total = fcount_total + fcount;
if (fcount ~= binary_length)
    error('Could not write CBF image-data.');
end
    
% write binary-end signature
[fcount] = fprintf(fid,'\r\n--CIF-BINARY-FORMAT-SECTION----\r\n;\r\n\r\n');
fcount_total = fcount_total + fcount;
if (fcount < 1)
    error('Could not write CBF binary-end signature');
end


% close output data file
fclose(fid);
if (debug_level >= 2)
    fprintf('%d data bytes written\n',fcount_total);
end

return;
