% Call function without arguments for instructions on how to use it

% Filename: $RCSfile: image_spec.m,v $
%
% $Revision: 1.19 $  $Date: 2012/08/23 15:14:57 $
% $Author:  $
% $Tag: $
%
% Description:
% display the current image data file
% (kind of online viewer)
%
% Note:
% information is passed on to this macro via a text file at a hard-wired
% location
%
% Dependencies:
% - image_show
% - common_header_value
%
%
% history:
%
% August 21th 2012
% Add a pause and retry for h5 files, the filename exists before the file
% is ready to read
%
% July 7th 2009:
% redraw figure if image size changes
%
% February 4th 2009:
% use valid-pixel mask for intensity sum calculation
%
% September 4th 2008:
% use new image_show functionality rather than displaying file date and
% exposure time from this macro
%
% August 28th 2008:
% pause one second in case of errors upon reading the filename from the
% spec exchange file
%
% June 24th 2008:
% update examples and include PCO CCD
%
% June 16th 2008:
% 1st version based on pilatus_spec

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

function [] = image_spec(det_no,varargin)
import beamline.pilatus_valid_pixel_roi
import io.common_header_value
import plotting.image_show
import plotting.image_spec
import utils.fopen_until_exists

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

% set default parameters
matlab_spec_filename = '~/Data10/specES1/internal/spec_matlab_print.dat';
% valid pixel mask
filename_valid_mask = '~/Data10/analysis/data/pilatus_valid_mask.mat';

% check number of input arguments
if ((nargin < 1) || (rem(no_of_in_arg,2) ~= 1))
    fprintf('Usage:\n')
    fprintf('%s(detector_number,[[<name>,<value>], ...]);\n',mfilename);
    fprintf('The optional name value pairs are:\n');
    fprintf('''MatlabSpecFilename'',<''dir and filename''>\n');
    fprintf('                                     default is %s\n',matlab_spec_filename);
    fprintf('''FilenameValidMask'',<''dir and filename ''>\n');
    fprintf('                                     default is %s\n',filename_valid_mask);
    fprintf('Additionally parameters of image_show are supported like:\n')
    fprintf('''RowFrom'',<0-max>                    region of interest definition, 0 for full frame\n');
    fprintf('''RowTo'',<0-max>                      region of interest definition, 0 for full frame\n');
    fprintf('''ColumnFrom'',<0-max>                 region of interest definition, 0 for full frame\n');
    fprintf('''ColumnTo'',<0-max>                   region of interest definition, 0 for full frame\n');
    fprintf('''AxisMin'',<value greater than 0>     intensity scaling\n');
    fprintf('''AxisMax'',<value greater than 0>     intensity scaling\n');
    fprintf('''ROI'',[ <ColumnFrom> <RowFrom> <ColumnTo> <RowTo> ]\n');
    fprintf('                                     region of interest definition of all four coordinates together\n');
    fprintf('Please call image_show without parameters for a complete list.\n')
    fprintf('The integral intensity is calculated for the specified region of interest and masks out\n');
    fprintf('invalid pixels via the valid pixel mask.\n');
    fprintf('\n');
    fprintf('\n');
    fprintf('Examples:\n');
    fprintf('%s;\n',mfilename);
    fprintf('%s(2);\n',mfilename);
    fprintf('%s(2,''RowFrom'',100,''RowTo'',800,''ColumnFrom'',200,''ColumnTo'',1200,''AxisMin'',1,''AxisMax'',1e6);\n',mfilename);
    fprintf('%s(2,''AxisMax'',1e4);\n',mfilename);
    fprintf('\n');
    fprintf('Recognized detector numbers are:\n');
    fprintf(' 1 - Pilatus 2M\n');
    fprintf(' 2 - Pilatus 300k\n');
    fprintf(' 3 - Pilatus 100k\n');
    fprintf(' 4 - Eiger single chip (256x256 pixels)\n');
    fprintf(' 5 - PX4 unit (Amptek)\n');
    fprintf(' 6 - MCS (SIS VME module)\n');
    fprintf(' 7 - FLI CCD\n');
    fprintf(' 8 - Mythen\n');
    fprintf(' 9 - Roper CCD\n');
    fprintf('10 - Andor CCD\n');
    fprintf('11 - PCO CCD\n');

    error('The optional parameters have to be specified as ''name'',''value'' pairs');
    return;
end

% for other than PILATUS detectors there is no default valid pixel mask
if (det_no > 2)
    filename_valid_mask = [];
end

% parse the variable input arguments not handled by image_show
vararg = cell(0,0);                    
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'MatlabSpecFilename' 
            matlab_spec_filename = value;
        case 'DetNo' 
            det_no = value;
        case 'FilenameValidMask' 
            filename_valid_mask = value;
        otherwise
            vararg{end+1} = name;
            vararg{end+1} = value;
    end
end

% add some parameters for image_show
vararg(13:(end+12)) = vararg;
vararg{ 1} = 'RetryReadSleep';
vararg{ 2} = 0.1;
vararg{ 3} = 'RetryReadMax';
vararg{ 4} = 10;
vararg{ 5} = 'MessageIfNotFound';
vararg{ 6} = 0;
vararg{ 7} = 'ErrorIfNotFound';
vararg{ 8} = 0;
vararg{ 9} = 'DisplayFilename';
vararg{10} = 0;
vararg{11} = 'ImageHandle';
image_handle_pos = 12;
vararg{image_handle_pos} = 0;

% set some default values for the plot window
% set(0, 'DefaultAxesfontsize', 12);
% set(0, 'DefaultAxeslinewidth', 1, 'DefaultAxesfontsize', 12);
% set(0, 'DefaultLinelinewidth', 1);

last_time_stamp = 0;
last_filename = ' ';
last_no_of_el_read = 0;
disp_ct = -1;
idle_ct = 0;
date_old = ' ';
date_new = ' ';

image_handle = 0;

plot_ind_in = 1;
plot_ind_out = 1;
plot_ind_max = 5;
plot_filename = cell(plot_ind_max,1);
for (ind = 1:plot_ind_max)
    plot_filename{ind} = '';
end

retry_read = 0;

last_date_ct = ' ';
last_frame_size = [];

% load the valid pixel mask ind_valid
if ((~isempty(filename_valid_mask)) && (exist(filename_valid_mask,'file')))
    fprintf('loading the valid pixel mask %s\n',filename_valid_mask);
    load(filename_valid_mask);
else
    fprintf('valid pixel mask not found, integral intensity includes hot pixels\n');
    valid_mask = [];
end


% periodically check the spec communication file
pause on;
while (true)
    err = 0;
    % open the plot-information file for read-only access
    [fid] = fopen_until_exists(matlab_spec_filename,'RetryReadSleep',10.0);
    
    % read up to including the line with the needed information
    for (ind=1:det_no)
        line = fgetl(fid);
        if (~ischar(line))
            fprintf('Could not read filename information for detector %d\n',det_no);
            err = 1;
            break;
        end
    end  
    if (err == 0)
        if (isempty(line))
            fprintf('Filename information for detector %d is empty\n',det_no);
            err = 1;
        end
    end

    % close the plot-information file
    if (fid > 0)
        if (fclose(fid) ~= 0)
            fprintf('Error upon closing %s\n',matlab_spec_filename);
            err = 1;
        end
    end
    
    % parse the line
    if (err == 0)
        det_name = [];
        time_stamp = [];
        filename = [];
        [det_name] = char(sscanf(line,'%[^:] %*f ''%*s'''));
        [time_stamp] = sscanf(line,'%*s %f ''%*s''');
        [filename] = char(sscanf(line,'%*s %*f ''%[^'']')');
        if ((isempty(det_name)) || (isempty(time_stamp)) || (isempty(filename)))
            fprintf('Could not parse filename information line for detector %d\n',det_no);
            err = 1;
        end
    else
        pause(1);
    end
    
    
    % store the filename, if it is new
    if ((err == 0) && (time_stamp ~= last_time_stamp))
        last_time_stamp = time_stamp;
        plot_filename{plot_ind_in} = filename;
        % suppress double names like multiple image_ct files
        if ((plot_ind_in == plot_ind_out) || (~strcmp(filename,last_filename)))
            plot_ind_in = plot_ind_in +1;
            if (plot_ind_in > plot_ind_max)
                plot_ind_in = 1;
            end
        end
        last_filename = filename;
    end

    % handle pending plots
    if (plot_ind_in ~= plot_ind_out)
        err = 0;

        if (retry_read == 0)
            fprintf('plotting %s   ',plot_filename{plot_ind_out});
        end

        % reuse the previous image handle for updating the figure
        vararg{image_handle_pos} = image_handle;
        
        if det_no == 4  % Add a default H5 data group for the Eiger
            auxno = numel(vararg);
            vararg{auxno+1} = 'H5Location';
            vararg{auxno+2} = '/eh5/images';
        end
        
        try
            [frame,image_handle] = image_show(plot_filename{plot_ind_out},vararg);
        catch
            fprintf('\nFailed reading, pausing 0.5 second and trying again\n')
            pause(0.5)
            [frame,image_handle] = image_show(plot_filename{plot_ind_out},vararg);
        end
        % redraw figure, if the figure size changed
        if ((~isempty(last_frame_size)) && ...
            (~isempty(frame)) && (~isempty(frame.data)) && ...
            (any(size(frame.data) ~= last_frame_size)))
            fprintf('image size changed, redrawing   ')
            set(gca, 'XLimMode','auto', 'YLimMode','auto');
            hold off;
            clf;
            image_handle = 0;
            vararg{image_handle_pos} = image_handle;
            [frame,image_handle] = image_show(plot_filename{plot_ind_out},vararg);
        end

        % do not change the zoom level upon updates
        set(gca, 'XLimMode','manual', 'YLimMode','manual');
        
        if (~isempty(frame.data))
            % set window title
            set(gcf,'Name','image_spec');
            % get and display date string from file
            date_new = common_header_value(frame.header{1},...
                frame.extension{1},'Date');
            % check, if an old file has been displayed
            is_ct = length(strfind(plot_filename{plot_ind_out},'_ct.'));
            if (((strcmp(date_new,date_old)) && ...
                 (last_no_of_el_read == frame.no_of_el_read{1})) || ...
                ((is_ct ~= 0) && (strcmp(date_new,last_date_ct))))
                % this is still the same image, retry
                retry_read = retry_read +1;
                pause(0.1);
            else
                % this is a new image, reset flags, increase counter
                date_old = date_new;
                last_no_of_el_read = frame.no_of_el_read{1};
                if (is_ct ~= 0)
                    last_date_ct = date_new;
                end
                retry_read = 0;
                plot_ind_out = plot_ind_out +1;
                if (plot_ind_out > plot_ind_max)
                    plot_ind_out = 1;
                end
   
                ind_valid = 1:numel(frame.data);
                if (isstruct(valid_mask))
                    % cut out the region of interest from the valid pixel mask
                    if ((isempty(last_frame_size)) || (sum(size(frame.data) ~= last_frame_size)) ~= 0)
                        last_frame_size = size(frame.data);
                        % cut out the valid mask for the currently active
                        % detector modules
                        valid_mask_cut = pilatus_valid_pixel_roi(valid_mask,...
                            'RoiSize',frame.img_full_size{1});
                        % from this cut out the region of interest, 
                        % if specified
                        valid_mask_cut = pilatus_valid_pixel_roi(valid_mask_cut,vararg);
                    end

                    ind_valid = valid_mask_cut.indices;
                end
                fprintf('integral %.3e counts\n',sum(frame.data(ind_valid)));
            end
        else
            % the image is not available, retry
            retry_read = retry_read +1;
            err = 1;
            pause(0.1);
        end
        
        if (retry_read > 10)
            fprintf('  could not load file\n');
            date_old = date_new;
            retry_read = 0;
            plot_ind_out = plot_ind_out +1;
            if (plot_ind_out > plot_ind_max)
                plot_ind_out = 1;
            end                
        end
        drawnow;
        idle_ct = 0;
    else
        idle_ct = idle_ct +1;
        if (idle_ct > 20) 
            drawnow;
            idle_ct = 0;
        end
    end
    
    pause(0.2);
end


