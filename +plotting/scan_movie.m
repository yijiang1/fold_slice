% Call function without arguments for instructions on how to use it

% Filename: $RCSfile: scan_movie.m,v $
%
% $Revision: 1.4 $  $Date: 2010/11/12 15:26:15 $
% $Author:  $
% $Tag: $
%
% Description:
% plot a sequence of 2D data frames, optionally save them as .avi file
%
% Note:
% Call without arguments for a brief help text.
%
% Dependencies: 
% - image_read
%
% history:
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

function [vararg_remain] = scan_movie(filenames,varargin)
import io.image_info
import io.image_read
import plotting.image_show
import utils.find_files


% set default values
% figure number for display
fig_no = 2;
% do not display time in title
display_time = 0;
% do not write AVI files
write_avi = 0;
% avi filename
avi_filename = '~/Data10/scan_movie.avi';
% AVI frames per second
avi_fps = 5;
% capture the full figure
capture_figure = 1;
% screen frames per second, 0 for no wait
screen_fps = 0;
% search masks rather than filenames are specified 
is_fmask = 1;
% use Linux/Unix find command to interprete filename mask
use_find = 1;

% check minimum number of input arguments
if (nargin < 1)
    fprintf('Usage:\n')
    fprintf('[vararg_remain]=%s(<filenames> [[,<name>,<value>] ...]);\n',...
            mfilename);
    fprintf('The optional <name>,<value> pairs are:\n');
    fprintf('''FigNo'',<integer value>              figure number for data display, default is %d\n',fig_no);
    fprintf('''WriteAVI'',<0-no, 1-yes>             write the image frames as AVI file, default is %d\n',write_avi);
    fprintf('''AVIfilename'',<filename>             default is %s\n',avi_filename);
    fprintf('''AVI_fps'',<value in Hz>              default is %.1f\n',avi_fps);
    fprintf('''CaptureFigure'',<0-axes,1-figure>    capture image or full figure, default is %d\n',capture_figure);
    fprintf('''Screen_fps'',<value in Hz>           pause if update is faster than this rate / 0 for no pause, default is %.1f\n',screen_fps);
    fprintf('''IsFmask'',<0-no,1-yes>               interprete the filename(s) as search mask that may include wildcards, default %d',is_fmask);
    fprintf('''UseFind'',<0-no, 1-yes>              use Linux/Unix command find to interprete the filename mask, default is %d\n',use_find);
    fprintf('Additional <name>,<value> pairs recognized by image_show can be specified. Please call image_show for an overview\n');
    fprintf('\n');
    fprintf('Examples:\n');
    fprintf('%s({''file1.cbf'',''file2.cbf''});\n',mfilename);
    fprintf('%s(''dir/*.cbf'');\n',mfilename);
    fprintf('%s(''dir/*.cbf'', ''WriteAVI'',1);\n',mfilename);
    fprintf('%s(''dir/*.cbf'',''ColumnFrom'',700,''ColumnTo'',1100,''RowFrom'',500,''RowTo'',1000,''Screen_fps'',0.5,''LogScale'',0,''AxisMax'',1e3);\n',mfilename);
    fprintf('%s(''dir/*.cbf'',''WriteAVI'',1,''AVIfilename'',''test.avi'');\n',mfilename);
    error('At least the filenames have to be specified as input parameter.');
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
vararg = cell(1,2);
% do not display current time by default, can be overruled with parameters
vararg{1} = 'DisplayTime';
vararg{2} = 0;
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'FigNo' 
            fig_no = value;
        case 'WriteAVI'
            write_avi = value;
        case 'AVIfilename'
            avi_filename = value;
        case 'AVI_fps' 
            avi_fps = value;
        case 'Screen_fps' 
            screen_fps = value;
        case 'IsFmask' 
            is_fmask = value;
        case 'UseFind' 
            use_find = value;
        case 'CaptureFigure'
            capture_figure = value;
        otherwise
            vararg{end+1} = name; %#ok<AGROW>
            vararg{end+1} = value; %#ok<AGROW>
    end
end

% pass some parameters to image_show
vararg{end+1} = 'FigNo';
vararg{end+1} = fig_no;
image_handle = 0;
vararg{end+1} = 'ImageHandle';
vararg{end+1} = image_handle;

% initialize the list of unhandled parameters
vararg_remain = cell(0,0);

% convert the filename to a cell array to use the same loop for single and
% multiple file names
if (~iscell(filenames))
    filenames = { filenames };
end

% convert single filename to cell array to ease handling
if (~is_fmask)
    fnames = struct('name',cell(1,1),'isdir',cell(1,1));
    fnames.isdir = false;
end

% loop over all filenames or filename masks
frame_no = 0;
tic;
previous_time = toc;
file_ind_max = length(filenames);
figure(fig_no);
hold off;
clf;
avi_obj = [];
for (file_ind = 1:file_ind_max)
    filename = filenames{file_ind};
    if (is_fmask)
        [data_dir,fnames,vararg_remain] = find_files( filename, 'UseFind',use_find );
        if (length(fnames) < 1)
            fprintf('No matching files found for %s.\n',filename);
            continue;
        end
    else
        fnames.name = filename;
    end
    for (fname_ind =1:length(fnames))
        if (~fnames(fname_ind).isdir)
            frame_no = frame_no +1;
            % initialize avi object in case this is the first frame
            if ((write_avi) && (isempty(avi_obj)))
                fprintf('Opening video output file %s, %.1f fps\n',...
                    avi_filename,avi_fps);
                avi_obj=VideoWriter(avi_filename);
                avi_obj.FrameRate=avi_fps;
                open(avi_obj);
            end
            % this does currently not work in case of several frames per
            % file
            vararg{end} = image_handle;
            filename = [data_dir fnames(fname_ind).name];
            fi = image_info(filename);
            for (frame_ind = 1:fi.no_of_frames)
                if (frame_ind == 1)
                    vararg_mod = vararg;
                    vararg_mod{length(vararg_mod)+1} = 'FrameNumber';
                    frame_ind_pos = length(vararg_mod)+1;
                    vararg_mod{frame_ind_pos} = frame_ind;
                    [frame,image_handle] = image_show(filename,vararg_mod);
                    % hand over the already read data to image_show for the
                    % display of the following frames
                    vararg_mod{length(vararg_mod)+1} = 'Frame';
                    vararg_mod{length(vararg_mod)+1} = frame;
                else
                    vararg_mod{frame_ind_pos} = frame_ind;
                    image_show(filename,vararg_mod);
                end
            
                drawnow;

                % add frame to avi file
                if (write_avi)
                    if (capture_figure)
                        avi_frame = getframe(fig_no);
                    else
                        avi_frame = getframe(gca);
                    end
                    writeVideo(avi_obj,avi_frame);
                end
                % pause in the unlikely case that the screen update is too fast
                if ((~isempty(screen_fps)) && (screen_fps > 0))
                    time_now = toc;
                    time_remain = 1.0/screen_fps - time_now + previous_time;
                    if (time_remain > 0)
                        pause(time_remain);
                        previous_time = 1.0/screen_fps + previous_time;
                    else
                        previous_time = time_now;
                    end
                end
            end
        end
    end
end

% close avi file, if necessary
if (write_avi)
    close(avi_obj); %#ok<NASGU>
end
