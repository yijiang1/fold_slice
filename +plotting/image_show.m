% Call function without arguments for instructions on how to use it

% Filename: $RCSfile: image_show.m,v $
%
% $Revision: 1.18 $  $Date: 2011/05/29 13:52:28 $
% $Author:  $
% $Tag: $
%
% Description:
%  read and display a 2D data frame
%
% Note:
% Call without arguments for a brief help text.
%
% Dependencies:
% - image_read
%
%
% history:
%
% April 16th 2009:
% bug-fix in histogram scaling, the intensity base level was not properly
% taken into account, 
% auto-scaling is activated by a separate parameter to allow for arbitrary
% scales
%
% February 23rd 2009:
% add histogram scaling
%
% September 4th 2009:
% add DisplayFtime and DisplayExptime parameters, 
% use new rowcol_from field of th eframe structure rather than handling the
% region of interest specification here,
% remove leading home directory and beamline specific path from displayed
% filename
%
% February 23rd 2009:
% add ImageHandle parameter for updating existing plots,
% earlier the XScange and YScange parameters have been added
%
% August 28th 2008:
% plot 1D data using plot
%
% May 9th 2008: 
% use image_read and named parameters, other changes since the last
% history entry, rename from pilatus_show to image_show
%
% November 2006: 1st version

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

function [frame,image_handle] = image_show(filename,varargin)
import io.*
import plotting.image_show
import utils.default_parameter_value


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

% set default values
fig_no = default_parameter_value(mfilename,'FigNo');
fig_clear = default_parameter_value(mfilename,'FigClear');
auto_scale = default_parameter_value(mfilename,'AutoScale');
axis_min = default_parameter_value(mfilename,'AxisMin');
axis_max = default_parameter_value(mfilename,'AxisMax');
hist_scale = default_parameter_value(mfilename,'HistScale');
log_scale = default_parameter_value(mfilename,'LogScale');
x_scale =default_parameter_value(mfilename,'XScale');
y_scale =default_parameter_value(mfilename,'YScale');
x_offs =default_parameter_value(mfilename,'XOffs');
display_colorbar = default_parameter_value(mfilename,'ColorBar');
display_axes = default_parameter_value(mfilename,'Axes');
display_time = default_parameter_value(mfilename,'DisplayTime');
display_ftime = default_parameter_value(mfilename,'DisplayFtime');
display_exptime = default_parameter_value(mfilename,'DisplayExptime');
color_map = default_parameter_value(mfilename,'ColorMap');
bgr_data = default_parameter_value(mfilename,'BgrData');
image_handle = default_parameter_value(mfilename,'ImageHandle');
frame_number = default_parameter_value(mfilename,'FrameNumber');
frame = [];

% check minimum number of input arguments and 
% check number of input arguments
if ((nargin < 1) || (rem(no_of_in_arg,2) ~= 1))
    image_read_help('ext',mfilename,'Examples',0,'ExtensionReturned',1);
    fprintf('''BgrData'',<data array>               subtract this background from the data read\n');
    fprintf('''FrameNumber'',<integer>              in case of multiple frames per file select one or 0 for the average of all frames, default is %d\n',frame_number);
    fprintf('''FigNo'',<figure number>              display figure in this window, default is %d\n',fig_no);
    fprintf('''FigClear'',<0-no, 1-yes>             clear figure before plotting, default is %d (only used if no image handle is specified)\n',fig_clear);
    fprintf('''ImageHandle'',<handle number>        update the specified image which may be faster and can be used to keep constant zoom level\n');
    fprintf('''AutoScale'',<[0-no/1-yes 0-no/1-yes]> intensity auto scaling of lower/upper limit, default is [%d %d]\n',auto_scale(1),auto_scale(2));
    fprintf('''AxisMin'',<value>                    intensity scaling if not in auto-scale mode, default is %.3e\n',axis_min);
    fprintf('''AxisMax'',<value>                    intensity scaling if not in auto-scale mode, default is %.3e\n',axis_max);
    fprintf('''AxisMinMax'',<[ min max]>            specify both min and max value\n');
    fprintf('''HistScale'',<[low high]>             in case of auto scaling scale the data to see the portion of the intensity values from low to high, default is [ %.2f %.2f]\n',hist_scale(1),hist_scale(2));
    fprintf('''LogScale'',<0-no,1-yes>              default is %d\n',log_scale);
    fprintf('''XScale'',<scaling factor>            default is %.3e\n',x_scale);
    fprintf('''YScale'',<scaling factor>            default is %.3e\n',y_scale);
    fprintf('''XOffs'',<scaling factor>             default is %.3e\n',x_offs);
    fprintf('''ColorBar'',<0-no,1-yes>              display colorbar, default is %d\n',display_colorbar);
    fprintf('''ColorMap'',<''map-name''>              choose colormap, default is ''%s''\n',color_map);
    fprintf('''Axes'',<0-no,1-yes>                  display axes, default is %d\n',display_axes);
    fprintf('''TitleString'',<''text''>               fixed part of the title\n');
    fprintf('''DisplayTime'',<0-no,1-yes>           default is yes\n');
    fprintf('\n');
    fprintf('\n');
    fprintf('Examples:\n');
    fprintf('[frame,image_handle]=%s(''~/Data10/pilatus/image_1_ct.cbf'');\n',mfilename);
    fprintf('[frame,image_handle]=%s(''~/Data10/pilatus/image_1_ct.cbf'');\n',mfilename);
    fprintf('[frame,image_handle]=%s(''~/Data10/pilatus/image_1_ct.cbf'',''RowFrom'',100,''RowTo'',800,''ColumnFrom'',200,''ColumnTo'',1200,''AxisMin'',1,''AxisMax'',1e6);\n',mfilename);
    fprintf('[frame,image_handle]=%s(''~/Data10/pilatus/image_1_ct.cbf'',''AxisMax'',1e4);\n',mfilename);
    fprintf('\n');
    fprintf('The returned structure has the fields data, header and extension.\n');
    if (nargin < 1)
        error('At least the filename has to be specified as input parameter.');
    else
        error('The optional parameters have to be specified as ''name'',''value'' pairs');
    end
end

% parse the variable input arguments not handled by or on purpose not
% passed to the image_read routine
vararg = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'FigNo' 
            fig_no = value;
        case 'FigClear'
            fig_clear = value;
        case 'ImageHandle'
            image_handle = value;
        case 'AxisMin' 
            axis_min = value;
        case 'AxisMax' 
            axis_max = value;
        case 'AxisMinMax'
            axis_min = value(1);
            axis_max = value(2);
        case 'AutoScale'
            auto_scale = value;
        case 'HistScale' 
            hist_scale = value;
        case 'LogScale' 
            log_scale = value;
        case 'XScale'
            x_scale = value;
        case 'YScale'
            y_scale = value;
        case 'XOffs'
            x_offs = value;
        case 'ColorBar'
            display_colorbar = value;
        case 'ColorMap'
            color_map = value;
        case 'Axes' 
            display_axes = value;
        case 'DisplayTime'
            display_time = value;
        case 'DisplayFtime'
            display_ftime = value;
        case 'DisplayExptime'
            display_exptime = value;
        case 'BgrData'
            bgr_data = value;
        case 'FrameNumber'
            frame_number = value;
        case 'Frame'
            frame = value;
        otherwise
            vararg{end+1} = name; %#ok<AGROW>
            vararg{end+1} = value; %#ok<AGROW>
    end
end

% get the home directory path
[stat,homedir] = system('echo ~');
if (stat ~= 0)
    homedir = '';
else
    if (~isempty(homedir))
        homedir = homedir(1:end-1);
    end
end

% remove a beamline specific part
[stat,username] = system('echo $USER');
if (stat == 0)
    if (~isempty(username))
        username = username(1:end-1);
    end
    std_path = [ '/sls/X12SA/Data10/' username ];
    if ((length(filename) > length(std_path)) && ...
        (strcmp(filename(1:length(std_path)),std_path)))
        filename = [ '~/Data10' filename(length(std_path)+1:end) ];
    end
end

% read the frame
if (isempty(frame))
    [frame] = image_read(filename,vararg);
end

if (isempty(frame.data))
     fprintf('%s: frame empty (file %s)\n',mfilename,filename);
     return;
end


if (axis_max <= axis_min)
    error('AxisMax must be greater than AxisMin');
end


% subtract background
if (~isempty(bgr_data))
    if (size(bgr_data) ~= size(frame.data))
        error('Background data dimension does not match data dimension.');
    end
    frame.data = frame.data - bgr_data;
end

% define axes range
x = x_scale * ...
    frame.rowcol_from{1}(2):(frame.rowcol_from{1}(2)+size(frame.data,2)-1) + ...
    x_offs;
y = y_scale * ...
    frame.rowcol_from{1}(1):(frame.rowcol_from{1}(1)+size(frame.data,1)-1);


% select the plot window
if ((isempty(image_handle)) || (image_handle <= 0))
    if (gcf ~= fig_no) 
      figure(fig_no);
    end
    if (fig_clear)
        hold off;
        clf;
    end
end


% for dat files (like MCS data) display the last data set as a 1D plot
if ((strcmp(frame.extension{1},'dat')) && (ndims(frame.data) > 1))
    x = x_scale * (1:size(frame.data,1)) + x_offs;
    y = y_scale;
    end_ind = floor(frame.no_of_el_read{1} / length(x)) * length(x);
    begin_ind = end_ind - length(x) +1;
    if (begin_ind < 1)
        begin_ind = 1;
    end
    frame_plot = frame.data(begin_ind:end_ind);
else
    if (ndims(frame.data) > 2)
        if (frame_number == 0)
            frame_plot = squeeze(mean(frame.data,3));
        else
            frame_plot = squeeze(frame.data(:,:,frame_number));
        end
    else
        frame_plot = frame.data;
    end
end


% display the frame
if (log_scale)
    frame_plot(frame_plot < 1e-15) = 1e-15;
end
ax_min = axis_min;
ax_max = axis_max;

if ((length(x) == 1) || (length(y) ==1))
    if (length(y) == 1)
        if (log_scale)
            semilogy(x,frame_plot);
        else
            plot(x,frame_plot);
        end

        xlabel('x [ pixel ]','FontSize',12);
        x_range = [ x(1) x(end) ];
    else
        if (log_scale)
            semilogy(y,frame_plot);
        else
            plot(y,frame_plot);
        end
        xlabel('y [ pixel ]','FontSize',12);
        x_range = [ y(1) y(end) ];
    end
    
    % axis scaling etc. 
    ylabel('intensity','FontSize',12);
    as = axis;
    as(1:2) = x_range;
    % auto scale lower/upper limit if specified
    if (auto_scale(1))
        ax_min = min(min(frame_plot));
    end
    if (auto_scale(2))
        ax_max = max(max(frame_plot));
    end
    as(3) = ax_min;
    as(4) = ax_max;
    if (as(4) <= as(3))
        as(4) = as(3) +1;
    end
    axis(as);
else
    if (log_scale)
        ax_min = log10( axis_min );
        ax_max = log10( axis_max );
        frame_plot = log10( double(frame_plot) );
    end

    if ((isempty(image_handle)) || (image_handle <= 0))
        % initialize a new plot
        image_handle = imagesc(x,y,frame_plot);
        axis xy;
        axis equal;
        axis tight;
        if (~display_axes)
            axis off;
        end
        jet_mod = jet;
        for (i = 1:8)
          jet_mod(i,3) = i/8;
        end
        if (isempty(color_map))
            colormap(jet_mod);
        else
            colormap(color_map);
        end
        if (display_colorbar)
            colorbar;
        end
        xlabel('x [ pixel ]','FontSize',12);
        ylabel('y [ pixel ]','FontSize',12);
        zlabel('intensity');
    else
        % update an existing plot
        set(image_handle,'CData', frame_plot );
    end
    
    % auto scale, if specified
    if (nnz(auto_scale))
        % use one million intensity bins over the intensity range
        int_min = min(min(frame_plot));
        int_max = max(max(frame_plot));
        int_step = (int_max-int_min) / 2^20;
        % calculate the bin indices from the intensities
        hist_ind = round((frame_plot(:)-int_min) / int_step) +1;
        % increase each bin an index is pointing to by one
        hist_bins = zeros(1,max(hist_ind)-min(hist_ind)+1);
        for (ind_hist = 1:numel(hist_ind))
            hist_bins(hist_ind(ind_hist)) = hist_bins(hist_ind(ind_hist)) +1;
        end
        % calculate the cumulative sum 
        hist_bins = cumsum(hist_bins);
        % find the indices to the bins with more pixels than the threshold
        n_low = find(hist_bins > hist_scale(1) * numel(frame_plot),1,'first');
        if (n_low >= length(hist_bins))
            n_low = length(hist_bins) -1;
        end
        n_high = find(hist_bins > hist_scale(2) * numel(frame_plot),1,'first');
        if (n_high <= n_low)
            n_high = n_low +1;
        end
        if (n_high > length(hist_bins))
            n_high = length(hist_bins);
        end
        
        % use the corresponding intensities to scale the image
        if (auto_scale(1))
            ax_min = (n_low-1) * int_step + int_min;
        end
        if (auto_scale(2))
            ax_max = (n_high-1) * int_step + int_min;
        end

%         for (hist_iter = 1:4)
%             hist_bins = int_min:((int_max-int_min)/1e3):int_max;
%             [no_of_int_values] = hist(frame.data(:),hist_bins);
%             % calculate the cumulative sum and find the position above the
%             % lower and upper number of pixels threshold
%             cumsum_no_of_int_values = cumsum(no_of_int_values);
%             n_low = find(cumsum_no_of_int_values > hist_scale(1) * numel(frame_plot),1,'first') -1;
%             if (n_low >= length(hist_bins))
%                 n_low = length(hist_bins) -1;
%             end
%             if (n_low < 1)
%                 n_low = 1;
%             end
%             n_high = find(cumsum_no_of_int_values > hist_scale(2) * numel(frame_plot),1,'first');
%             if (n_high <= n_low)
%                 n_high = n_low +1;
%             end
%             if (n_high > length(hist_bins))
%                 n_high = length(hist_bins);
%             end
% 
%             % sufficiently fine resolution is reached
%             if ((n_high > 5) && (n_high < length(hist_bins)))
% %                 fprintf('%d: %d %d',hist_iter,n_low,n_high);
%                 break;
%             end
%             
%             % change to a finer histogram spacing
%             int_min = hist_bins(n_low);
%             int_max = hist_bins(n_high);
%         end

%         % use the corresponding intensities to scale the image
%         if (ax_min <= 0)
%             ax_min = hist_bins(n_low);
%         end
%         if (ax_max <= 0)
%             ax_max = hist_bins(n_high);
%         end

%         frame_med = medfilt2(frame_plot,[ 3 3 ]);
%         if (ax_min <= 0)
%             ax_min = 0.1 * mean(mean(frame_med));
%         end
%         if (ax_max <= 0)
%             ax_max = max(max(frame_med));
%         end
    end
    caxis( [ ax_min ax_max ] );

end


% compile title string:

% replace home directory by tilde
if ((length(filename) > length(homedir)) && ...
    (strcmp(filename(1:length(homedir)),homedir)))
    title_str = [ '~' filename(length(homedir)+1:end) ];
else
    title_str = filename;
end
% escape underscore characters from LaTeX-like use
title_str = strrep(title_str,'_','\_');
if (log_scale)
    title_str = [ title_str ' (log.)' ];
end
% add current time
if (display_time)
    clock_now = clock();
    title_str = [ title_str '   ' ...
                  num2str(clock_now(4)) ':' num2str(clock_now(5),'%02.0f') ...
                  ':' num2str(clock_now(6),'%02.0f') ];
end
% get and display date string from file
if (display_ftime) 
    title_str = [ title_str '\newlinefile date: ' ...
        common_header_value(frame.header{1},frame.extension{1},'Date') ];
end

% get and display exposure time
if (display_exptime)
    title_str = [ title_str ', exposure time = ' ...
        num2str(common_header_value(frame.header{1},frame.extension{1},...
                'ExposureTime'),'%.3f') ' sec' ];
end

% add the frame number to the title string
if (frame_number == 0)
    title_str = [ title_str ', frame average' ];
else
    title_str = [ title_str sprintf(', frame %d',frame_number) ];
end

% display the compiled title string
title(title_str,'FontSize',14);


return;
