% Call function without arguments for instructions on how to use it

% Filename: $RCSfile: speread.m,v $
%
% $Revision: 1.1 $  $Date: 2008/06/10 17:05:14 $
% $Author:  $
% $Tag: $
%
% Description:
% Macro for reading SPE CCD data files written by WinView from 
% Roper Scientific / PI Acton
%
% Note:
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
% May 9th 2008: 1st version
%

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

function [frame,vararg_remain] = speread(filename,varargin)
import io.*
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
    image_read_sub_help(mfilename,'SPE');
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
max_header_length = 4100;

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

% some sanity checks
if (fdat(4099:4100) ~= [ 85; 85 ])
    error('%s: no header end signature found',filename);
    return;
end


% return some selected information from the header as lines of a cell array
% and extract some variables needed for the data extraction
% version 2.5 is assumed, otherwise some values may be wrong
frame.header = cell(0,0);
xdim = double(typecast(fdat(43:44),'uint16'));
frame.header{end+1} = sprintf('xdim %.0f',xdim);

ydim = double(typecast(fdat(657:658),'uint16'));
frame.header{end+1} = sprintf('ydim %.0f',ydim);

exposure_sec = double(typecast(fdat(11:14),'single'));
frame.header{end+1} = sprintf('exposure %.3f',exposure_sec);

frame.header{end+1} = sprintf('xDimDet %d',typecast(fdat(7:8),'uint16'));
frame.header{end+1} = sprintf('yDimDet %d',typecast(fdat(19:20),'uint16'));

frame.header{end+1} = sprintf('date %s',char(fdat(21:30)'));
frame.header{end+1} = sprintf('ExperimentTimeLocal %s:%s:%s',...
    char(fdat(173:174)'),char(fdat(175:176)'),char(fdat(177:178)'));

data_type = fdat(109);
frame.header{end+1} = sprintf('datatype %d',data_type);
switch (data_type)
    case 0
        bytes_per_pixel = 4;
        cast_type = 'single';
    case 1
        bytes_per_pixel = 4;
        cast_type = 'int32';
    case 2
        bytes_per_pixel = 2;
        cast_type = 'int16';
    case 3
        bytes_per_pixel = 2;
        cast_type = 'uint16';
    otherwise
        error('unknown data type %d',data_type);
end

frame.header{end+1} = sprintf('BackGrndApplied %d',typecast(fdat(151:152),'uint16'));

frame.header{end+1} = sprintf('gain %d',typecast(fdat(199:200),'uint16'));

% frame.header{end+1} = sprintf('comments %s',char(fdat(201:600)'));

% number of accumulations
frame.header{end+1} = sprintf('lavgexp %ld',typecast(fdat(669:672),'uint32'));

frame.header{end+1} = sprintf('ReadoutTime %f',typecast(fdat(673:676),'single'));

frame.header{end+1} = sprintf('sw_version %s',char(fdat(689:703)'));

frame.header{end+1} = sprintf('type %d',typecast(fdat(705:706),'uint16'));

frame.header{end+1} = sprintf('flatFieldApplied %d',typecast(fdat(707:708),'uint16'));

num_frames = double(typecast(fdat(1447:1450),'uint32'));
frame.header{end+1} = sprintf('NumFrames %.0f',num_frames);

file_header_ver = double(typecast(fdat(1993:1996),'single'));
frame.header{end+1} = sprintf('file_header_ver %.2f',file_header_ver);


% extract the data frames
data_bytes = bytes_per_pixel * xdim * ydim * num_frames;
if (fcount-4100 ~= data_bytes)
    error('%d bytes expected, %d found',data_bytes,fcount-4100);
end
frame.data = reshape(typecast(fdat(4101:end),cast_type),...
    xdim, ydim, num_frames );

% conversion to standard view on PI-SCX4300 data at the SLS/cSAXS beamline
% if (~original_orientation)
%     frame.data = permute(frame.data,[2 1 3]);
% end

% return as double precision floating point
frame.data = double(frame.data);

 