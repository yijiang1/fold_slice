% Call function without arguments for a detailed explanation of its use

% Filename: $RCSfile: fliread.m,v $
%
% $Revision: 1.4 $  $Date: 2008/10/03 13:50:04 $
% $Author:  $
% $Tag: $
%
% Description:
% read a data file in the format stored by the program ccdfli.c
%
% Note:
% The image files have the extension raw and the file format is home
% defined. 
% Call without arguments for a brief help text.
%
% Dependencies:
% - image_read_set_default
% - fopen_until_exists
% - get_hdr_val
%
%
% history:
%
% Pctober 3rd 2008: 
% add new end-of-header signature search for version 1.20 files
%
% October 1st 2008: Exchange width and height in reshape command
%
% May 9th 2008: adapt to call from image_read
%
% November 3, 2005: include new fields of data format 1.1: 
% exposure time Spec, exposure time measured, monitor counts
%
% October 2005: include optional from-to line reading
%
% March 2005: 1st version

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

function [frame,vararg_remain] = fliread(filename, varargin)
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
    image_read_sub_help(mfilename,'raw');
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
max_header_length = 1024;

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

% return the complete header as lines of a cell array
frame.header = char_to_cellstr( char(fdat(1:max_header_length)'),1 );

version_no = get_hdr_val(frame.header,'% version','%f',1);
if ((version_no ~= 1.00) && (version_no ~= 1.10) && (version_no ~= 1.20))
     fprintf('%s: File version number %.2f may not be supported\n',...
	     mfilename,version_no);
end

[frameHeight,line_number] = get_hdr_val(frame.header,'% rows','%d',1);
[frameWidth,line_number]  = get_hdr_val(frame.header,'% columns','%d',1);

if (version_no < 1.20)
    if (version_no == 1.10)
        [monCounts,line_number]  = get_hdr_val(frame.header,'% monitorcounts','%d',1);
    end

    % cut off non-header lines
    frame.header = frame.header(1:line_number);

    % find start of data
    eol_ind = regexp(char(fdat(1:max_header_length)'),'\n');
    data_start = eol_ind(line_number) +1;
else
    eoh_signature = sprintf('%% EOH%c%c',10,26);
    end_of_header_pos = ...
        strfind( fdat(1:min(max_header_length,length(fdat)))',...
             eoh_signature );
    data_start = end_of_header_pos + length(eoh_signature);
end


% calculate end of data (should be end of file)
data_end = data_start + frameWidth * frameHeight *2 -1;
if (data_end > fcount)
    error('%d bytes read but %d are needed',fcount,data_end);
end
if (data_end ~= fcount)
     fprintf('%s warning: %d bytes read vs. %d needed\n',mfilename,...
             fcount,data_end);
end
% cut out frame data
frame.data = double( reshape(typecast(fdat(data_start:data_end),'uint16'), ...
    frameWidth,frameHeight) );

% conversion to standard view on FLI-CCD data at the SLS/cSAXS beamline
% (to be determined)
% if (~original_orientation)
% %     frame = flipud(frame');
%     frame.data = frame.data';
% end

% add the file modification date to the header
dir_entry = dir(filename);
frame.header{end+1} = [ 'FileTimestamp ' dir_entry.date ];
