% Call function without arguments for instructions on how to use it

% Filename: $RCSfile: image_read.m,v $
%
% $Revision: 1.17 $  $Date: 2013/01/25 10:23:23 $
% $Author: $
% $Tag: $
%
% Description:
% Macro for reading image data formats used at the SLS / cSAXS beamline. 
% The data are returned in double precision floating point format. 
%
% Note:
% Call without arguments for a brief help text.
%
% Dependencies: 
% - edfread
% - cbfread
% - hdf5read
% - fliread
% - speread
% - char_to_cellstr
%
% history:
%
% September 30th 2010:
% add a call to hdf5read_main
%
% June 5th 2009:
% disable UhandledParError before calling sub-macros
%
% January 16th 2009:
% adapt to image_orient returning the complete structure rather than just
% the data array
%
% November 19th 2008:
% add reading of Matlab files
%
% September 5th 2008:
% skip further processing for a frame if it was not possible to read it
%
% September 4th 2008:
% add the rowcol-from field to the frames structure as origin information
%
% June 19th 2008: adapt call to image_orient
% 
% May 16th 2008: send variable arguments through find files
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
%   using the “cSAXS matlab package” developed by the CXS group
%   and Computing Department, Paul Scherrer Institut, Switzerland.” 
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

function [frames,vararg_remain] = image_read(filenames,varargin)
import io.*
import io.HDF.*
import io.CBF.*
import plotting.*
import utils.char_to_cellstr
import utils.default_parameter_value
import utils.find_files
import utils.fopen_until_exists

% initialize return arguments
frames = struct('data',[], ...
    'img_full_size',[], 'rowcol_from',[], ...
    'no_of_el_read', [], ...
    'header',[], 'filename',[], 'extension', []);

% check minimum number of input arguments
if (nargin < 1)
    image_read_help('ext',mfilename);
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

% convert the filename to a cell array to use the same loop for single and
% multiple file names
if (~iscell(filenames))
    filenames = { filenames };
end

% hdf5 files are read via a separate sub-routine
% if length(filenames) == 1 % accept only one file name
    filename = filenames{1};
    [~, ~, ext] = fileparts(filename);
    
    if any(strcmp(ext, {'.h5', '.hdf5', '.nxs', '.cxs'}))
        if length(filenames) == 1
            frames = hdf5read(filename, varargin);
        else
            frames = hdf5read(filenames, varargin);
        end
        vararg_remain = [];
        
        return % image_read ends here for hdf5 image files
    end
% end

% set default values for the variable input arguments:
% default data type for the returned frames
data_type = default_parameter_value(mfilename,'DataType');
% recognize file type by file name extension
force_file_type = default_parameter_value(mfilename,'ForceFileType');
% from/to row 0 means all rows
row_from = default_parameter_value(mfilename,'RowFrom');
row_to = default_parameter_value(mfilename,'RowTo');
% from/to column 0 means all lines
column_from = default_parameter_value(mfilename,'ColumnFrom');
column_to = default_parameter_value(mfilename,'ColumnTo');
% determine default orientation based on the file name extension
orient_by_extension = default_parameter_value(mfilename,'OrientByExtension');
% filename is actually a mask that may include wildcards
filename_is_fmask = default_parameter_value(mfilename,'IsFmask');
% display file name of the file to be loaded
display_filename = default_parameter_value(mfilename,'DisplayFilename');
% variable to load from Matlab files
matlab_var = default_parameter_value(mfilename,'MatlabVar');

% exit with an error message if unhandled named parameters are left at the
% end of this macro
unhandled_par_error = 1;

% parse the variable input arguments
vararg = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'DataType'
            if (~ischar(value))
                error('The DataType must be string defining a valid Matlab data type.');
            end
            data_type = value;
        case 'ForceFileType'
            force_file_type = lower(value);
        case 'MatlabVar' 
            matlab_var = value;
        case 'RowFrom' 
            row_from = value;
        case 'ROI'
            if (length(value) ~= 4)
                error('The ROI parameter needs a vector of length four as argument.');
            end
            column_from = value(1);
            row_from = value(2);
            column_to = value(3);
            row_to = value(4);
        case 'RowTo' 
            row_to = value;
        case 'ColumnFrom' 
            column_from = value;
        case 'ColumnTo' 
            column_to = value;
        case 'OrientByExtension' 
            orient_by_extension = value;
        case 'UnhandledParError'
            unhandled_par_error = value;
        case 'IsFmask' 
            filename_is_fmask = value;
        case 'DisplayFilename' 
            display_filename = value;
        otherwise
            vararg{end+1} = name; %#ok<AGROW>
            vararg{end+1} = value; %#ok<AGROW>
    end
end


% initialize the list of unhandled parameters
vararg_remain = cell(0,0);

% loop over all specified file names
file_ind_max = length(filenames);
store_ind = 1;
for (file_ind=1:file_ind_max)
    filename = filenames{file_ind};
    vararg_remain = vararg;

    % in case of file name mask get a list of all matching file names
    data_dir = '';
    if (filename_is_fmask)
        % sub macros must not complain about unknown arguments
        vararg_remain{end+1} = 'UnhandledParError';
        vararg_remain{end+1} = 0;
        
        [data_dir, fnames, vararg_remain] = ...
            find_files( filename, vararg_remain );
    else
        fnames = struct('name',filename);
    end
    
    for (sub_file_ind = 1:length(fnames))
        % pick out the current filename
        filename = [ data_dir fnames(sub_file_ind).name ];
        
        % check for minimum filename length
        if (length(filename) < 5)
            error([ mfilename ': invalid filename ' filename ]);
        end

        if (isempty(force_file_type))
            % get the extension from the last three to four characters
            extension = lower(filename((end-4):end));
            pos = strfind(extension,'.');
            if (length(pos) < 1)
                error([ mfilename ': invalid extension in ' filename ]);
            end
            extension = extension(pos(end)+1:end);
        else
            % the file name extension is ignored since the file type is
            % forced to a specific one
            extension = force_file_type;
        end

        if (display_filename)
            fprintf('loading %s\n',filename);
        end
        
        if ((strcmp(extension,'dat')) || ...
            (strcmp(extension,'tif')) || (strcmp(extension,'tiff')) || ...
            (strcmp(extension,'mat')))
            % open the file to support functionality like 
            % wait-until-exists
            [fid,vararg_remain] = fopen_until_exists(filename,vararg);
            if (fid >= 0)
              fclose(fid);
            end        
        end
        
        % interprete file in the format indicated by the filename extension
        switch extension
            case 'cbf'
                [frame,vararg_remain] = cbfread(filename,vararg_remain);
            case {'hdf5', 'h5', 'nxs', 'cxs'}
                [frame,vararg_remain] = hdf5read_main(filename,vararg_remain);
            case 'dat' 
                [frame,vararg_remain] = datread(filename,vararg_remain);
            case 'edf'
                [frame,vararg_remain] = edfread(filename,vararg_remain);
            case {'mar', 'mccd'}
                [frame,vararg_remain] = marread(filename,vararg_remain);
            case 'mat' 
                tmp_data = load(filename);
                frame.data = tmp_data.(matlab_var);
                frame.header = {};
                % No header information are available. 
                % Fake exposure time information to avoid problems in other
                % macros. 
                frame.header{end+1} = 'Exposure_time 1.0';
                % add the file modification date to the header
                dir_entry = dir(filename);
                frame.header{end+1} = [ 'DateTime ' dir_entry.date ];
            case 'raw'
                [frame,vararg_remain] = fliread(filename,vararg_remain);
            case 'spe'
                [frame,vararg_remain] = speread(filename,vararg_remain);
            case {'tif', 'tiff'}
                % reading higher bit depths than 16bit needs a sufficiently
                % up-to-date Matlab version
                frame.data = imread(filename,'tif');
                hdr = imfinfo(filename);
                frame.header = {};
                if (isfield(hdr,'ImageDescription'))
                    frame.header = char_to_cellstr(hdr.ImageDescription);
                end
                if ((isfield(hdr,'Model')) && ...
                    (strcmp(hdr.Model(1:7),'PILATUS')))
                    frame.header{end+1} = [ 'DateTime ' hdr.DateTime ];
                    frame.header{end+1} = [ 'Software ' hdr.Software ];
                    frame.header{end+1} = [ 'Model ' hdr.Model ];
                else
                    % the exposure time is not available
                    if (isfield(hdr,'exptimesec'))
                        frame.header{end+1} = [ 'Exposure_time' hdr.exptimesec ];
                    else
                        frame.header{end+1} = 'Exposure_time 1.0';
                    end
                    if (isfield(hdr,'DateTime'))
                        frame.header{end+1} = [ 'DateTime ' hdr.DateTime ];
                    else
                        if (isfield(hdr,'FileModDate'))
                            frame.header{end+1} = [ 'DateTime ' hdr.FileModDate ];
                        else
                            % add the file modification date to the header
                            dir_entry = dir(filename);
                            frame.header{end+1} = ...
                              [ 'DateTime ' dir_entry.date ];
                        end
                    end
                end
            otherwise
                error([ 'unknown extension of ' filename ]);
        end

        % the remaining code is not needed if no file was read
        if (isempty(frame.data))
            continue;
        end
        % determine the orientation from the filename extension
        vararg_remain_prev = vararg_remain;
        vararg_remain = cell(1,length(vararg_remain)+2);
        vararg_remain(3:end) = vararg_remain_prev;
        vararg_remain{1} = 'OrientByExtension';
        vararg_remain{2} = orient_by_extension;
        % set the extension since image_orient is called prior
        % to defining the return variables frames 
        frame.extension = cell(1,1);
        frame.extension{1} = extension;
        
        % orient image
        [frame,vararg_remain] = image_orient(frame,vararg_remain);
        
        % cut out region of interest
        full_size = size(frame.data);
        if ((row_from > 0) || (row_to > 0) ||... 
            (column_from > 0) || (column_to > 0))
            if (row_from <= 0)
                row_from = 1;
            end
            if (row_from > size(frame.data,1))
                error('The RowFrom specification is beyond the maximum value of %d',...
                    size(frame.data,1));
            end
            if (row_to <= row_from)
                row_to = full_size(1);
            end
            if (row_to > full_size(1))
                error('The RowTo specification is beyond the maximum value of %d',...
                    full_size(1));
            end

            if (column_from <= 0)
                column_from = 1;
            end
            if (column_from > full_size(2))
                error('The ColumnFrom specification is beyond the maximum value of %d',...
                    full_size(2));
            end
            if (column_to <= column_from)
                column_to = full_size(2);
            end
            if (column_to > size(frame.data,2))
                error('The ColumnTo specification is beyond the maximum value of %d',...
                    full_size(2));
            end
            frame.data = frame.data(row_from:row_to,column_from:column_to,:);
        end

        % initialize the return array with the now known dimensions
        if (store_ind == 1)
            % in case of file name masks or multiple images in one data
            % file the final array dimensions can only be estimated
            init_guess = file_ind_max -1 + length(fnames);
            frames.data = zeros( [ size(frame.data) init_guess ], data_type );
            frames.rowcol_from = cell(1,init_guess);
            frames.no_of_el_read = cell(1,init_guess);
            frames.img_full_size = cell(1,init_guess);
            frames.filename = cell(1,init_guess);
            frames.extension = cell(1,init_guess);
            frames.header = cell(1,init_guess);
        end


        % store the frame(s)
        if ((size(frame.data,1) ~= size(frames.data,1)) || ...
            (size(frame.data,2) ~= size(frames.data,2)))
            error('Expected frame dimension is %d x %d, frame read has %d x %d',...
                size(frames.data,2),size(frames.data,1),...
                size(frame.data,2),size(frame.data,1));
        end

        % in some cases multiple frames are stored in a single file
        if (ndims(frame.data) == 4)
            store_ind_to = (store_ind+size(frame.data,4)-1);
            frames.data(:,:,:,store_ind:store_ind_to) = cast(frame.data,data_type);
        else
            store_ind_to = (store_ind+size(frame.data,3)-1);
            frames.data(:,:,store_ind:store_ind_to) = cast(frame.data,data_type);
        end
        
        % store the filename, extension and the header in the return argument
        for (ind=store_ind:store_ind_to)
            frames.img_full_size{ind} = full_size;
            frames.rowcol_from{ind} = [ row_from column_from ];
            if (isfield(frame,'no_of_el_read'))
                frames.no_of_el_read{ind} = frame.no_of_el_read;
            else
                frames.no_of_el_read{ind} = size(frame.data,3);
            end
            % zero means no ROI, i.e., starting at point (1,1)
            if (frames.rowcol_from{ind}(1) < 1)
                frames.rowcol_from{ind}(1) = 1;
            end
            if (frames.rowcol_from{ind}(2) < 1)
                frames.rowcol_from{ind}(2) = 1;
            end
            frames.filename{ind} = filename;
            frames.extension{ind} = extension;
            frames.header{ind} = frame.header;
        end
        
        % update index to free space in the output arrays
        store_ind = store_ind_to + 1;

        % exit in case of unhandled named parameters, if this has not been switched
        % off
        if ((unhandled_par_error) && (~isempty(vararg_remain)))
            vararg_remain
            error('Not all named parameters have been handled.');
        end
        % restore parameters for next iteration
        vararg_remain = vararg_remain_prev;        
    end
end

% resize the output arrays in case the initial size is too large
store_ind = store_ind -1;
if (size(frames.data,3) > store_ind)
    frames.data = frames.data(:,:,1:store_ind);
    frames.img_full_size = frames.img_full_size(1:store_ind);
    frames.rowcol_from = frames.rowcol_from(1:store_ind);
    frames.filename = frames.filename(1:store_ind);
    frames.extension = frames.extension(1:store_ind);
    frames.header = frames.header(1:store_ind);
end

function frames = hdf5read(filename, params)
import utils.fopen_until_exists
import io.image_default_orientation
import io.HDF.hdf5_load
import io.image_orient
import utils.find_files

frames = struct('data', [], 'img_full_size', [], 'rowcol_from', [], ...
    'no_of_el_read', [], 'header',[], 'filename',[], 'extension', []);



p = inputParser;
p.KeepUnmatched = true;
p.FunctionName = 'image_read';

addParameter(p, 'H5Location', '/');
addParameter(p, 'ReadAttr', false);

addParameter(p, 'FrameRange', [1, Inf], @(x) isvector(x) && numel(x) <= 2 && isnumeric(x));
addParameter(p, 'RowFrom', 1, @(x) isscalar(x) && isnumeric(x));
addParameter(p, 'RowTo', Inf, @(x) isscalar(x) && isnumeric(x));
addParameter(p, 'ColumnFrom', 1, @(x) isscalar(x) && isnumeric(x));
addParameter(p, 'ColumnTo', Inf, @(x) isscalar(x) && isnumeric(x));

addParameter(p, 'OrientByExtension', 1, @(x) isscalar(x) && (isnumeric(x) || islogical(x)));
addParameter(p, 'Orientation', [0, 0, 0], @(x) isvector(x) && numel(x) == 3 && isnumeric(x));
addParameter(p, 'InvertOrientation', 0, @(x) isscalar(x) && isnumeric(x));
addParameter(p, 'Transpose', 0, @(x) isscalar(x) && isnumeric(x));
addParameter(p, 'FlipLR', 0, @(x) isscalar(x) && isnumeric(x));
addParameter(p, 'FlipUD', 0, @(x) isscalar(x) && isnumeric(x));
addParameter(p, 'CatDim', -1, @(x) isscalar(x) && isnumeric(x));

addParameter(p, 'DisplayFilename', 1, @(x) isscalar(x) && isnumeric(x));

% No support for filename wildcards, i.e. ('IsMask', 1)
addParameter(p, 'IsFmask', true, ...
    @islogical);%@(x) assert(~x, 'Filename wildcards, i.e. (''IsFmask'', 1), are not supported for hdf5 files'));

parse(p, params{:});
r = p.Results;
vararg_remain = [fieldnames(p.Unmatched)'; struct2cell(p.Unmatched)'];

if (r.IsFmask)    
    [data_dir, fnames, vararg_remain] = ...
        find_files( filename, vararg_remain );
    filename = [];
    for ii=1:length(fnames)
        filename{ii} = fullfile(data_dir, fnames(ii).name);
    end
end


if iscell(filename)
    dir_entry = dir(filename{1});
    frames.header{1} = {['DateTime ' dir_entry.date], 'Exposure_time 0.0'};
    [~, ~, frames.extension{1}] = fileparts(filename{1});
    frames.extension{1} = frames.extension{1}(2:end); % truncate a leading dot
    frames.filename{1} = '*multiple_frames*';
else
    dir_entry = dir(filename);
    frames.header{1} = {['DateTime ' dir_entry.date], 'Exposure_time 0.0'};
    [~, frames.filename{1}, frames.extension{1}] = fileparts(filename);
    frames.extension{1} = frames.extension{1}(2:end); % truncate a leading dot
end


inputs_orient = {'OrientByExtension', r.OrientByExtension, 'InvertOrientation', r.InvertOrientation};

if r.OrientByExtension
    orient_vec = image_default_orientation(frames.header{1}, frames.extension{1});
    do_transpose = orient_vec(1);
    do_fliplr = orient_vec(2);
    do_flipud = orient_vec(3);
    
    % Warn user if they also provided either of 'Orientation', 'Transpose', 'FlipLR' or 'FlipUD'
    if ~all(ismember({'Orientation', 'Transpose', 'FlipLR', 'FlipUD'}, p.UsingDefaults))
        warning(['Default hdf5 image orientation is potentially modified by either ' ...
            '''Orientation'', or any of ''Transpose'', ''FlipLR'', or ''FlipUD'' parameters. ' ...
            'To supress this warning set ''OrientByExtension'' to 0.']);
    end
    
else
    do_transpose = 0;
    do_fliplr = 0;
    do_flipud = 0;
end

if ~ismember({'Orientation'}, p.UsingDefaults)
    do_transpose = r.Orientation(1);
    do_fliplr = r.Orientation(2);
    do_flipud = r.Orientation(3);
    
    inputs_orient = [inputs_orient, {'Orientation', r.Orientation}];
    
    % Warn user if they also provided either of 'Transpose', 'FlipLR' or 'FlipUD'
    if ~all(ismember({'Transpose', 'FlipLR', 'FlipUD'}, p.UsingDefaults))
        warning(['Image orientation specified via ''Orientation'' parameter is potentially ' ...
            'modified by either ''Transpose'', ''FlipLR'', and/or ''FlipUD''. ' ...
            'To supress this warning use either ''Orientation'' or a combination of ' ...
            '''Transpose'', ''FlipLR'', and/or ''FlipUD'' parameters.']);
    end
end

if r.Transpose || r.FlipLR || r.FlipUD
    do_transpose = r.Transpose;
    do_fliplr = r.FlipLR;
    do_flipud = r.FlipUD;
    inputs_orient = [inputs_orient, {'Transpose', r.Transpose}];
    inputs_orient = [inputs_orient, {'FlipLR', r.FlipLR}];
    inputs_orient = [inputs_orient, {'FlipUD', r.FlipUD}];
end

inputs = {};
if ~isempty(r.H5Location) && ischar(r.H5Location)
    inputs{end+1} = r.H5Location;
end

if r.ReadAttr
    inputs{end+1} = '-sa';
end

% Add slicing indexes if a user specified any of them
if ~all(ismember({'FrameRange', 'RowFrom', 'RowTo', 'ColumnFrom', 'ColumnTo'}, ...
        p.UsingDefaults))
    
    % Support 0's as start/end index -> full left/right range
    if numel(r.FrameRange) == 1
        r.FrameRange(2) = r.FrameRange(1);
    end
    if r.FrameRange(1) == 0; r.FrameRange(1) = 1; end
    if r.FrameRange(2) == 0; r.FrameRange(2) = Inf; end
    if r.RowFrom == 0; r.RowFrom = 1; end
    if r.RowTo == 0; r.RowTo = Inf; end
    if r.ColumnFrom == 0; r.ColumnFrom = 1; end
    if r.ColumnTo == 0; r.ColumnTo = Inf; end
    
    % Adjust range values according to the consequent image orientation procedure
    if r.InvertOrientation % Transpose -> FlipLR/FlipUD
        if do_transpose; [r.ColumnTo, r.ColumnFrom, r.RowTo, r.RowFrom] = ...
                deal(r.RowTo, r.RowFrom, r.ColumnTo, r.ColumnFrom); end
        if do_fliplr; [r.ColumnTo, r.ColumnFrom] = deal(-r.ColumnFrom, -r.ColumnTo); end
        if do_flipud; [r.RowTo, r.RowFrom] = deal(-r.RowFrom, -r.RowTo); end
        
    else % FlipLR/FlipUD -> Transpose
        if do_fliplr; [r.ColumnTo, r.ColumnFrom] = deal(-r.ColumnFrom, -r.ColumnTo); end
        if do_flipud; [r.RowTo, r.RowFrom] = deal(-r.RowFrom, -r.RowTo); end
        if do_transpose; [r.ColumnTo, r.ColumnFrom, r.RowTo, r.RowFrom] = ...
                deal(r.RowTo, r.RowFrom, r.ColumnTo, r.ColumnFrom); end
    end
    
    % Form the input
    inputs{end+1} = {[r.RowFrom, r.RowTo],[r.ColumnFrom, r.ColumnTo], r.FrameRange};
end
if iscell(filename)
    
    if (r.DisplayFilename)
        fprintf('loading %s\n', filename{1});
    end
    tmp = frames;
    tmp.data = hdf5_load(filename{1}, inputs{:});
    % Orient image frame(s)
    [tmp, vararg_remain] = image_orient(tmp, [inputs_orient, vararg_remain]);
    
    % let's handle the 2D (or 3d with singleton) case first
    if isnumeric(tmp.data) && ndims(tmp.data==3) && size(tmp.data,3)==1
        frames.data = zeros([size(tmp.data(:,:,1)) length(filename)*size(tmp.data,3)]);
        frames.data(:,:,1:size(tmp.data,3)) = tmp.data;
        if length(filename)>1
            for frame=2:length(filename)
                if (r.DisplayFilename)
                    fprintf('loading %s\n', filename{frame});
                end
                tmp.data = hdf5_load(filename{frame}, inputs{:});
                
                % Orient image frame(s)
                [tmp, vararg_remain] = image_orient(tmp, [inputs_orient, vararg_remain]);
                frames.data(:,:,frame) = tmp.data;
            end
        end
    else
        % in case of more than 2 dimensions, concatenate along the
        % specified dimension, or return a cell array
        if length(filename)>1
            frames.data{1} = tmp.data;
            framedim = ndims(frames.data{1});
            for frameID=2:length(filename)
                if (r.DisplayFilename)
                    fprintf('loading %s\n', filename{frameID});
                end
                frames.data{frameID} = hdf5_load(filename{frameID}, inputs{:});
                if ndims(frames.data{frameID})~=framedim
                    framedim = -1;
                end
                % Orient image frame(s)
                [frames, vararg_remain] = image_orient(frames, [inputs_orient, vararg_remain]);
            end
            if framedim > 0
                try
                    if r.CatDim == -1
                        frames.data = cat(ndims(frames.data{1}),frames.data{:});
                    else
                        frames.data = cat(r.CatDim, frames.data{:});
                    end
                catch
                    warning('Failed to concatenate frames.')
                end
            end
                    
        else
            frames.data = tmp.data;
        end
    end
else
    % Support wait-until-exist functionality
    [fid, vararg_remain] = fopen_until_exists(filename, vararg_remain(:));
    if fid >= 0
        fclose(fid);
    else
        % Silently exit if a file was not found and ('ErrorIfNotFound', 0)
        return;
    end
    
    
    if (r.DisplayFilename)
        fprintf('loading %s\n', filename);
    end
    frames.data = hdf5_load(filename, inputs{:});
    
    % Orient image frame(s)
    [frames, vararg_remain] = image_orient(frames, [inputs_orient, vararg_remain]);
    
end

% If its a dataset then fill these fields for further showing with
% image_show.m or image_spec.m
if isnumeric(frames.data)
    frames.img_full_size = {[size(frames.data, 1), size(frames.data, 2)]};
    frames.rowcol_from = {[r.RowFrom, r.ColumnFrom]};
    frames.no_of_el_read = {size(frames.data, 3)};
elseif (do_fliplr||do_flipud||do_transpose)
    warning(['H5Location points to a group, not a dataset. Orientation/OrientByExtension/Transpose/FlipUD/FlipLR will be ignored. \n '...
    'To remove this warning, set ''OrientByExtension'' to 0 and ''Orientation'' to [0 0 0]'])
end

% Show a warning message for unsupported input parameters
unmatched = vararg_remain(1:2:end);
if ~isempty(unmatched)
    warning('These input parameters are not supported for hdf5 files and will be ignored: %s', ...
        strjoin(unmatched, ', '));
end
