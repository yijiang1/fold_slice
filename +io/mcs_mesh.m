% Call function without arguments for instructions on how to use it

% Filename: $RCSfile: mcs_mesh.m,v $
%
% $Revision: 1.7 $  $Date: 2016/08/03 08:38:32 $
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

function [mcs_data, data_adjusted, pos_data] = mcs_mesh(first_scan_no,no_of_intervals,varargin)
import io.*
import utils.adjust_projection
import utils.fopen_until_exists
import utils.get_hdr_val

% initialize return arguments
mcs_data = [];
% legacy_sgalil = false;  % A flag that keeps track of different commands needed if the older legacy file of reading positions is used - 2019.04

% set default parameter
% plot this MCS channel
ch_to_plot = 4;
snake_scan=0;
fast_axis_x = 1;
% create the plot in this figure
fig_no = 123;
% exit with an error message if unhandled named parameters are left at the
% end of this macro
unhandled_par_error = 1;
% file name base
fname_base = '';
% first part of directory path
dir_base = '~/Data10/mcs/'; 
% scaling factors for the axes
x_scale = 1.0;
y_scale = 1.0;
%
axis_minmax = [];
% save resulting figure
figure_dir = '~/Data10/analysis/online/stxm/figures/';      
% save the resulting data
data_dir =  '~/Data10/analysis/online/stxm/data/';      
pos_file = '~/Data10/sgalil/S%05d.dat';
positions_only = false;

% check minimum number of input arguments
if (nargin < 2)
    fprintf('[mcs_data data_adjusted pos_data]=%s(<first scan no.>, <no. of line intervals> [[,<name>,<value>] ...]);\n',...
        mfilename);
    fprintf('The optional <name>,<value> pairs are:\n');
    fprintf('''ChToPlot'',<channel no.>              if greater than zero than this channel is plotted, default is %d\n',ch_to_plot);
    fprintf('''SnakeScan'',<0-no, 1-yes>             scan mode is a snake pattern, default is 0\n');
    fprintf('''FastAxisX'',<0-no, 1-yes>             fast scan axis is along x, default is 1\n');
    fprintf('''FigNo'',<figure number>               plot the data in this figure, 0 for no figure, default is %d\n',fig_no);
    fprintf('''XScale'',<scale factor>               scale the x-axis with this factor\n');
    fprintf('''YScale'',<scale factor>               scale the y-axis with this factor\n');
    fprintf('''AxisMinMax'',<[ min max]>             specify both min and max value\n');    
    fprintf('''DirBase'',<first part of data path>   including an ending slash, default is ''%s''\n',dir_base);
    fprintf('''FnameBase'',<first part of file name> in case of an empty string the current Unix user name is used followed by an underscore, default if ''%s''\n',fname_base);
    fprintf('''FigureDir'',''directory''               save the resulting plot in eps, jpeg and Matlab fig format, '''' for no saving, default is %s\n',figure_dir);
    fprintf('''DataDir'',''directory''                 save the resulting data as Matlab file, '''' for no saving, default is %s\n',data_dir);
    fprintf('''Pos_file'',''path and filename structure''   checks the flipping of images using the positions found in these files. It gives a warning if it detects the wrong fast axis and it flips the lines for snake scans.\n');
    fprintf('                                           If the filename structure includes an %% it tries to read positions assuming one dat file for each scan. If it does not include an %% then it uses spec_read.\n');
    fprintf('                                           Default is %s, set to empty =[ ] to avoid warnings \n',pos_file);
    fprintf('''Positions_only'',<0-no, 1-yes>        if you only use the routine to get and check positions of a scan but you don''t have a MCS scalar measurement, default is %d.\n', positions_only);
    fprintf('The data are in the format ''fast to slow axis'', i.e., the first dimension is the MCS channel,\n');
    fprintf('the second dimension is the exposure index, the third dimension is the fast axis of a mesh scan\n');
    fprintf('and the fourth dimension is the slow axis of a mesh scan.\n');
    fprintf('An optional second output ''data_adjusted'' can be obtained. This is a structure that has been adjusted and flipped according to the scan parameters in order to reflect\n');
    fprintf('  the sample physical orientation. One of the fields is ''transm'' which is the sample transmissivity. Other fields are positions_out, scan_num, and scan_point.\n');
    error('At least the number of the first scan and the number of line intervals need to be specified as input parameter.');
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
            no_of_in_arg = 1 + length(varargin);
        end
    end
end

% check number of input arguments
if (rem(no_of_in_arg,2) ~= 0)
    error('The optional parameters have to be specified as ''name'',''value'' pairs');
end

% parse the variable input arguments
vararg = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'ChToPlot'
            ch_to_plot = value;
        case 'SnakeScan'
            snake_scan = value;
        case 'FastAxisX'
            fast_axis_x = value;
        case 'FigNo'
            fig_no = value;
        case 'XScale'
            x_scale = value;
        case 'YScale'
            y_scale = value;
        case 'AxisMinMax'
            axis_minmax = value;
        case 'UnhandledParError'
            unhandled_par_error = value;
        case 'FnameBase'
            fname_base = value;
        case 'DirBase'
            dir_base = value;
        case 'DataDir'
            data_dir = value;
        case 'FigureDir'
            figure_dir = value;
        case 'Pos_file'
            pos_file = value;
        case 'Positions_only'
            positions_only = value;
        otherwise
            vararg{end+1} = name; %#ok<AGROW>
            vararg{end+1} = value; %#ok<AGROW>
    end
end


% initialize the list of unhandled parameters
vararg_remain = cell(0,0);


% get the current user name

    if (length(fname_base) < 1)
        [stat,usr]=unix('echo $USER');
        fname_base = [ sscanf(usr,'%s') '_' ];
    end


% initialize the output figure
if ~positions_only
    figure(fig_no);
    hold off;
    clf;
end

last_scan_no = first_scan_no + no_of_intervals;
store_ind = 1;
last_draw_time = clock;

% over all the scan lines
for scan_no = first_scan_no:last_scan_no
    dir = [ dir_base 'S' num2str(floor(scan_no/1000)*1000,'%05d') '-' ...
        num2str(floor(scan_no/1000)*1000+999,'%05d') '/S' num2str(scan_no,'%05d') '/' ];
    filename =  [ fname_base num2str(scan_no,'%05d') '.dat' ];
    
    % read the frame until all data are available
    ind_rep = 0;
    ind_max = 3;
    last_no_of_el_read = 0;
    if ~positions_only
        while (ind_rep < ind_max)
            frame = image_read([dir filename ], 'RetryReadSleep',10, ...
                'RetryReadMax',1, 'RetrySleepWhenFound',10);
            if ~isempty(frame.data)&&(frame.no_of_el_read{1} >= numel(frame.data))
                % the complete data set has been read
                break;
            end
            if ~isempty(frame.data)&&(frame.no_of_el_read{1} <= last_no_of_el_read)
                % no progress, increase timeout counter
                ind_rep = ind_rep +1;
            else
                ind_rep = 0;
            end
            if ~isempty(frame.no_of_el_read)
                last_no_of_el_read = frame.no_of_el_read{1};
            end
            if ~isempty(frame.header)
                exp_time = get_hdr_val(frame.header{1},'Exposure_time','%f',1);
            end
            wait_time = numel(frame.no_of_el_read)*exp_time;
            if (wait_time > 2.0)
                fprintf('%d/%d: frame incomplete (%d/%d), waiting %.1fs and retrying\n',...
                    ind_rep+1,ind_max,frame.no_of_el_read{1},numel(frame.data));
            end
            pause(wait_time);
        end
    end
    if ~isempty(pos_file)
        if contains(pos_file,'%')
            positions_from_spec = false;
            filepos = sprintf(pos_file,scan_no);
            try
                positions = beamline.read_position_file(filepos);
                positions.data(1,:) = positions.Avg_x;
                positions.data(2,:) = positions.Avg_y;
                numpts = numel(positions.Avg_x);
            catch
                warning('The reading of positions for sgalil did not work, now trying legacy mode in older sgalil position format')
                positions = image_read(filepos,'RetryReadSleep',10,'RetryReadMax',0);
                numpts = size(positions.data,2);
%                 legacy_sgalil = true;
            end
        else % If it does not contain % delimiter then we assume is a spec file
            positions_from_spec = true;
            if scan_no == last_scan_no % Only read spec positions once all scans are done, otherwise is too slow to read each time
                try
                    fprintf('Reading positions from spec file %s \n',pos_file)
                    positions = io.read_scan_positions_spec(pos_file,first_scan_no:last_scan_no,{'samx','samy'});
                catch
                    warning('io.read_scan_positions_spec failed, pausing 5 seconds and retrying')
                    pause(5);
                    positions = io.read_scan_positions_spec(pos_file,first_scan_no:last_scan_no,{'samx','samy'});
                end
            else
                positions.data = [];
            end
        end
    end
    % initialize return array
    if (scan_no == first_scan_no)
        if ~positions_only
            mcs_data = zeros(size(frame.data,1),size(frame.data,2),...
                size(frame.data,3) + 1 ,no_of_intervals+1);
            numpts = size(frame.data,3) + 1;
        else
            mcs_data = [];
            numpts = size(positions.data,2);
        end
        
        [scan_num, scan_point] = meshgrid(first_scan_no:last_scan_no,0:numpts-1);
        % Check positions of stage for flipping
        if ~isempty(pos_file)
            pos_data =  zeros(numpts, no_of_intervals+1, 2);
        else
            pos_data = [];
        end
    end
    
    if ~positions_only
        store_ind_to = store_ind + size(frame.data,4) - 1;
    else
        store_ind_to = store_ind;
    end
    %%% Special about MCS is that it does not take the first image because
    %%% it is triggered in a strange way, so here to match other
    %%% detectors we make the matrix one element larger and we replicate
    %%% the first value
    if ~positions_only
        mcs_data(:,:,2:end,store_ind:(store_ind+size(frame.data,4)-1)) = frame.data;
        mcs_data(:,:,1,store_ind:(store_ind+size(frame.data,4)-1)) = frame.data(:,1,2,:);
    end

    if ~isempty(pos_file)
        if ~positions_from_spec
            pos_data(:,store_ind,:) = positions.data.';
        else % if positions are from spec they are read only at the end and arranged in a structure so they need special handling
            if scan_no == last_scan_no
                pos_data =  zeros(size(positions(1).data,1),no_of_intervals+1,2);
                for ii = numel(positions)
                    pos_data(:,ii,:) = positions(ii).data;
                end
            end
        end
    end
    store_ind = store_ind_to +1;
    
    if ((scan_no == first_scan_no) || (scan_no == last_scan_no) ||  ...
        (etime(clock,last_draw_time) > 10))
        if (size(mcs_data,2) == 1)
            data_plot = squeeze(mcs_data(ch_to_plot,1,:,:));
        elseif ~isempty(mcs_data)
            data_plot = squeeze(mcs_data(ch_to_plot,:,:,1));
        else
            data_plot = [];
        end
        if (~positions_only)&&((size(data_plot,1) > 1) && (size(data_plot,2) > 1))
            %CHANGE
            [data_plot, positions_out] = adjust_projection(data_plot, snake_scan, fast_axis_x, pos_data);
            % 2D plot
            x_values = (1:size(data_plot,2)) * x_scale;
            y_values = (1:size(data_plot,1)) * y_scale;

            if (~isempty(axis_minmax))
                caxis(axis_minmax);
            end
            %plot the image
            figure(fig_no)
            imagesc(data_plot);
            axis xy;
            axis equal;
            axis tight;
            colormap gray;
            colorbar;
            title( [ fname_base ': #'  num2str(first_scan_no) ' -' num2str(scan_no)] );
            drawnow;
        elseif (~positions_only)
            % 1D plot if only one line has been read
            x_values = (1:length(data_plot)) * x_scale;
            plot(x_values,data_plot);
        else
            [data_plot, positions_out] = adjust_projection(data_plot, snake_scan, fast_axis_x, pos_data);
        end
        last_draw_time = clock;
    end
end

[scan_num]   = adjust_projection(scan_num,   snake_scan, fast_axis_x, pos_data);
[scan_point] = adjust_projection(scan_point, snake_scan, fast_axis_x, pos_data);

data_adjusted.transm = data_plot;
data_adjusted.positions_out = positions_out;
data_adjusted.scan_num = scan_num;
data_adjusted.scan_point = scan_point;

% file name for saving
filename = sprintf('stxm_scans_%05d-%05d_mcs',first_scan_no,...
    last_scan_no);

% save figures
if (~isempty(figure_dir))&&(~positions_only)
    figure(fig_no);

    % create output directories and write the plot in different formats
    if (~exist(figure_dir,'dir'))
        mkdir(figure_dir)
    end
    if ((figure_dir(end) ~= '/') && (figure_dir(end) ~= '\')) 
        figure_dir = [ figure_dir '/' ];
    end
    fprintf('output directory for figures is %s\n',figure_dir);
    
    subdir = [ figure_dir 'jpg/' ];
    if (~exist(subdir,'dir'))
        mkdir(subdir);
    end
    fprintf('saving %s.jpg\n',filename);
    print('-djpeg','-r300',[subdir filename '.jpg'] );
    
    subdir = [ figure_dir 'eps/' ];
    if (~exist(subdir,'dir'))
        mkdir(subdir);
    end
    fprintf('saving %s.eps\n',filename);
    print('-depsc','-r1200',[subdir filename '.eps'] );
    
    subdir = [ figure_dir 'fig/' ];
    if (~exist(subdir,'dir'))
        mkdir(subdir);
    end
    fprintf('saving %s.fig\n',filename);
    hgsave([subdir filename '.fig']);    
end

% save resulting data
if (~isempty(data_dir))
    if ((data_dir(end) ~= '/') && (data_dir(end) ~= '\')) 
        data_dir = [ data_dir '/' ];
    end
    
    % create output directory
    if (~exist(data_dir,'dir'))
        mkdir(data_dir)
    end
    
    % save data
    fprintf('saving %s.mat\n',[data_dir filename]);
    save([data_dir filename],'mcs_data','first_scan_no','last_scan_no','pos_data','data_adjusted');
end

return

