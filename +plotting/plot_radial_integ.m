% Call function without arguments for instructions on how to use it

% Filename: $RCSfile: plot_radial_integ.m,v $
%
% $Revision: 1.15 $  $Date: 2016/01/21 14:50:38 $
% $Author:  $
% $Tag: $
%
% Description:
% plot radially integrated intensities
%
% Note:
% Call without arguments for a brief help text.
% The integrated intensities should be calculated first using
% radial_integ.m.
%
% Dependencies:
% none
%
% history:
%
% November 25th 2011:
% bug-fix in the background subtraction, only subtract for available
% intensities, i.e., intensities not flagged as -1
%
% May 27th 2011:
% return all curves rather than just the last one plotted and allow to
% suppress plotting completely using FigNo 0 (to average data using this
% function)
%
% April 28th 2010:
% add plot as a function of angle option,
% use default_parameter_value
%
% February 19th 2009:
% average only positive intensities (i.e., valid pixels)
%
% September 4th 2009: add Axis parameter
%
% June 9th 2008: 1st documented version

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
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

function [x_values_returned,y_values_returned] = plot_radial_integ(filename_masks, varargin)
import utils.default_parameter_value
import utils.find_files
import utils.pixel_to_q

% set default values
fig_no = default_parameter_value(mfilename,'FigNo');
new_fig = default_parameter_value(mfilename,'NewFig');
clear_fig = default_parameter_value(mfilename,'ClearFig');
axis_scale = default_parameter_value(mfilename,'Axis');
sleep_time = default_parameter_value(mfilename,'SleepTime');
xlog = default_parameter_value(mfilename,'XLog');
ylog = default_parameter_value(mfilename,'YLog');
plot_q = default_parameter_value(mfilename,'PlotQ');
plot_angle = default_parameter_value(mfilename,'PlotAngle');
radius_range = default_parameter_value(mfilename,'RadiusRange');
filename_integ_masks = default_parameter_value(mfilename,'FilenameIntegMasks');
pixel_size_mm = default_parameter_value(mfilename,'PixelSize_mm');
det_dist_mm = default_parameter_value(mfilename,'DetDist_mm');
E_keV = default_parameter_value(mfilename,'E_keV');
inverse_nm = default_parameter_value(mfilename,'Inverse_nm');
q_mul_pow = default_parameter_value(mfilename,'QMulPow');
seg_avg = default_parameter_value(mfilename,'SegAvg');
seg_range = default_parameter_value(mfilename,'SegRange');
legend_mul_seg = default_parameter_value(mfilename,'LegendMulSeg');
point_avg = default_parameter_value(mfilename,'PointAvg');
point_range = default_parameter_value(mfilename,'PointRange');
bgr_filename = default_parameter_value(mfilename,'BgrFilename');
bgr_scale = default_parameter_value(mfilename,'BgrScale');
bgr_point = default_parameter_value(mfilename,'BgrPoint');
show_fig = 1;

% check minimum number of input arguments
if (nargin < 1)
    fprintf('\nUsage:\n');
    fprintf('[x,y]=%s(filename_mask,  [[,<name>,<value>] ...]);\n',mfilename);
    fprintf('x and y are optional outputs to get the data which is plotted\n')
    fprintf('filename_mask can be something like ''*_integ.mat'' or ''*_integ.txt'' or\n');
    fprintf('a cell array of filenames or filename masks like {''dir1/*.mat'',''dir2/*.mat''}.\n');
    fprintf('The optional <name>,<value> pairs are:\n');
    fprintf('''FigNo'',<figure number>              number of the figure for plotting the integrated intensities, default is %d\n',fig_no);
    fprintf('''NewFig'',<0-no, 1-yes>               open a new figure for each file, default is %d\n',new_fig);
    fprintf('''ClearFig'',<0-no, 1-yes for the first point,2-yes, always>\n');
    fprintf('                                     clear the figure before plotting, default is %d\n',clear_fig);
    fprintf('''Axis,<[ x_from x_to y_from y_to ]>   fixed scale for the plot\n')
    fprintf('''SleepTime'',<seconds>                wait time after each plot, default is %.3f\n',sleep_time);
    fprintf('''XLog'',<0-no, 1-yes>                 logarithmic scaling of the x-axis, default is %d\n',xlog);
    fprintf('''YLog'',<0-no, 1-yes>                 logarithmic scaling of the y-axis, default is %d\n',ylog);
    fprintf('''PlotQ'',<0-no, 1-yes>                plot as a function of momentum transfer q rather than pixel no., default is %d\n',plot_q);
    fprintf('''PlotAngle'',<0-no, 1-yes>            plot as a function of the azimuthal angle rather than q or the radius, default is %d\n',plot_angle);
    fprintf('''RadiusRange'',<vector or []>         for azimuthal plots the intensity over this radius range is averaged, default is [] for all radii\n');
    fprintf('''FilenameIntegMasks'',<filename>      Matlab file containing the integration masks, needed for normalization in case of averaging over radii, default is ''%s''\n',filename_integ_masks);
    fprintf('''QMulPow'',<value,or []>              multiply intensity with q to the power of this value, default is [ ] for no multiplication\n');
    fprintf('''Inverse_nm'',<0-no, 1-yes>           plot q in inverse nm rather than inverse Angstroem, default is %d\n',inverse_nm);
    fprintf('''PixelSize_mm'',<value in mm>         pixel size for q-calculation, default is %.3f mm. If the file contains a q-vector this parameter is ignored\n',pixel_size_mm);
    fprintf('''DetDist_mm'',<value in mm>           sample to detector distance for q-calculation, default is %.3f mm. If the file contains a q-vector this parameter is ignored\n',det_dist_mm);
    fprintf('''E_keV'',<value in keV>               x-ray energy for q-calculation, default is %.3f mm. If the file contains a q-vector this parameter is ignored\n',E_keV);
    fprintf('''SegAvg'',<0-no, 1-yes>               average over angular segments rather than plotting them with different line colours, default is %d\n',seg_avg);
    fprintf('''SegRange'',<vector or []>            segment range to plot, default is [] for all segments\n');
    fprintf('''LegendMulSeg'',<0-no, 1-yes>         show a legend in case of multiple segments being plotted, default is %d\n',legend_mul_seg);
    fprintf('''PointAvg'',<0-no,1-yes>              plot the average of all intensity curves in the file, which typically means the average of a scan line, default is %d\n',point_avg);
    fprintf('''PointRange'',<vector or []>          point range to plot, default is [] for all points in a file\n');
    fprintf('''BgrFilename'',<''filename''>           background to subtract from each intensity profile, must have the same dimensions the data have\n');
    fprintf('''BgrScale'',<value>                   scaling factor to apply to the backgroubnd data, default is %.3e\n',bgr_scale);
    fprintf('''BgrPoint'',<integer>                 point to use from the file BgrFilename, default is %d, use [] to subtract 1:1\n',bgr_point);
    fprintf('''ShowFigure'',<0-no,1-yes>              show the figure,  default is %d\n',show_fig);
    
    fprintf('Examples:\n');
    fprintf('%s(''~/Data10/analysis/integ/data1_integ.mat'');\n',mfilename);
    error('At least the filename mask has to be specified as input argument.');
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
% vararg = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'FigNo'
            fig_no = value;
        case 'NewFig'
            new_fig = value;
        case 'ClearFig'
            clear_fig = value;
        case 'Axis'
            if (isempty(value))
                continue;
            end
            if (length(value) ~= 4)
                error('Axis needs a vector with four components as argument.');
            end
            axis_scale = value;
        case 'SleepTime'
            sleep_time = value;
        case 'XLog'
            xlog = value;
        case 'YLog'
            ylog = value;
        case 'PlotQ'
            plot_q = value;
        case 'PlotAngle'
            plot_angle = value;
        case 'RadiusRange'
            radius_range = value;
        case 'FilenameIntegMasks'
            filename_integ_masks = value;
        case 'QMulPow'
            q_mul_pow = value;
        case 'Inverse_nm'
            inverse_nm = value;
        case 'PixelSize_mm'
            pixel_size_mm = value;
        case 'DetDist_mm'
            det_dist_mm = value;
        case 'E_keV'
            E_keV = value;
        case 'SegRange'
            seg_range = value;
        case {'SegSum', 'SegAvg' }
            seg_avg = value;
        case 'LegendMulSeg'
            legend_mul_seg = value;
        case {'PointSum', 'PointAvg'}
            point_avg = value;
        case 'PointRange'
            point_range = value;
        case 'BgrFilename'
            bgr_filename = value;
        case 'BgrScale'
            bgr_scale = value;
        case 'BgrPoint'
            bgr_point = value;
        case 'ShowFigure'
            show_fig = value;
        otherwise
            error('Do not know how to handle the parameter %s',name);
            %             vararg{end+1} = name;
            %             vararg{end+1} = value;
    end
end

% automatically disable averaging over angles in case of angular plots
if (plot_angle)
    if (seg_avg ~= 0)
        seg_avg = 0;
        fprintf('Disabling the averaging over azimuthal segments since a plot as the function of angle has been requested.\n');
    end
end

% The integration masks are needed for normalization in case of averaging
% over radii.
if ((plot_angle) && (exist(filename_integ_masks,'file')))
    fprintf('Loading the integration masks from %s\n',filename_integ_masks);
    integ_data = load(filename_integ_masks);
else
    integ_data = [];
end

% determine the 1D plot function to use
if ((~xlog) && (~ylog))
    plot_function = @plot;
end
if ((~xlog) && (ylog))
    plot_function = @semilogy;
end
if ((xlog) && (~ylog))
    plot_function = @semilogx;
end
if ((xlog) && (ylog))
    plot_function = @loglog;
end


% load the background data
if (~isempty(bgr_filename))
    fprintf('reading background data from %s\n',bgr_filename);
    [d_bgr] = load(bgr_filename);
end

% loop over all filename masks

if (isstruct(filename_masks))
    % allow for other macros handing over directly the integrated data
    % rather than a filename
    data_provided = 1;
    ind_mask_max = 1;
    d = filename_masks;
else
    % ease handling by ensuring that filename_masks is a cell array
    if (~iscell(filename_masks))
        filename_masks = { filename_masks };
    end
    
    ind_mask_max = length(filename_masks);
    data_provided = 0;
end

y_values_returned = [];
x_values_returned = [];

for (ind_mask = 1:ind_mask_max) %#ok<*NO4LP>
    
    if (data_provided)
        file_ind_max = 1;
    else
        filename_mask = filename_masks{ind_mask};
        
        %     % get data directory
        %     [data_dir] = fileparts(filename_mask);
        %     if ((~isempty(data_dir)) && (data_dir(end) ~= '/'))
        %         data_dir = [ data_dir '/' ];
        %     end
        
        % search matching filenames
        [ data_dir, fnames, vararg_remain ] = find_files( filename_mask );
        if (length(fnames) < 1)
            fprintf('No matching files found for %s.\n',filename_mask);
            continue;
        end
        
        % loop over all matching files
        file_ind_max = length(fnames);
    end
    
    if (new_fig)
        legend_str = zeros(file_ind_max,12);
    end
    
    first_plot = 1;
    
    for (file_ind=1:file_ind_max)
        if (show_fig)
        if (((new_fig) || (first_plot)) && (fig_no > 0))
            figure(fig_no);
            if (clear_fig)
                hold off;
                clf;
            else
                if (((strcmp(get(gca,'XScale'),'linear')) && (xlog)) || ...
                        ((strcmp(get(gca,'YScale'),'linear')) && (ylog)))
                    hold off;
                else
                    hold all;
                end
            end
        end
        if (new_fig)
            first_plot = 1;
        end
        end
        if (~data_provided)
            % skip sub directories
            if (fnames(file_ind).isdir)
                continue;
            end
            
            % read one data file
            filename = [ data_dir fnames(file_ind).name ];
            fprintf('reading %4d / %4d: %s\n',file_ind,file_ind_max,...
                filename);
            [d] = load(filename);
        else
            filename = 'online plot';
        end
        
        if (isa(d,'struct'))
            % in Matlab files a structure with the data is stored
            if (isfield(d,'radius'))
                % current name
                radius = d.radius;
            elseif (isfield(d,'r'))
                % old name
                radius = d.r;
            else
                radius = d.q;
            end
            I_all = d.I_all;
            
            if (~isempty(bgr_filename))
                % subtract background
                if (isempty(bgr_point))
                    % subtract for the available intensities (I >= 0) 1:1,
                    % e.g., a line from a line
                    ind_pos = intersect(find(I_all >= 0), find(d_bgr.I_all >= 0));
                    I_all = I_all(ind_pos) - bgr_scale * d_bgr.I_all(ind_pos);
                else
                    % determine the over the specified point range
                    % averaged background intensity,
                    % keep unavailable intensities flagged as -1
                    I_all_bgr = calc_I_point_avg(d_bgr.I_all,bgr_point);
                    % subtract the average background from all intensity
                    % distributions
                    for (ind_point = 1:size(I_all,3))
                        % only subtract available intensities, not the with
                        % -1 flagged invalid ones
                        I_point = I_all(:,:,ind_point);
                        ind_pos = intersect(find(I_point >= 0), find(I_all_bgr >= 0));
                        I_point(ind_pos) = I_point(ind_pos) - bgr_scale * I_all_bgr(ind_pos);
                        I_all(:,:,ind_point) = I_point;
                    end
                end
            end
        else
            % in text files the first column contains the radius, the rest
            % are intensities for the different segments
            if ((isa(d,'double')) && (size(d,2) >= 2))
                radius = d(:,1);
                I_all = d(:,2:end);
                if (~isempty(bgr_filename))
                    I_all = I_all - bgr_scale * d_bgr(:,2:end);
                end
            else
                error('Unknown data format.');
            end
        end
        
        % check if segments have been loaded and get the number of segments
        % available or specified via a command line parameter
        no_of_radii = size(I_all,1);
        no_of_segments = size(I_all,2);
        no_of_points = size(I_all,3);
        if (no_of_segments < 1)
            error('could not load %s',filename);
        end
        if ((plot_angle) && (no_of_segments < 2))
            error('At least two segments need to be present for an angular plot.');
        end
        % average over segments if specified
        if (seg_avg)
            no_of_segments = 1;
            if (isempty(seg_range))
                seg_range_use = 1:size(I_all,2);
            else
                seg_range_use = seg_range;
            end
            % average over all pixels with positive intensities, i.e., skip
            % negative intensities
            I_all_prev = I_all;
            I_all = zeros(no_of_radii,1,no_of_points);
            for (ind1=1:no_of_radii)
                for (ind3=1:no_of_points)
                    no_of_el = 0;
                    for (ind2=1:length(seg_range_use))
                        if (I_all_prev(ind1,seg_range_use(ind2),ind3) >= 0)
                            I_all(ind1,1,ind3) = I_all(ind1,1,ind3) + I_all_prev(ind1,seg_range_use(ind2),ind3);
                            no_of_el = no_of_el +1;
                        end
                    end
                    if (no_of_el > 1)
                        I_all(ind1,1,ind3) = I_all(ind1,1,ind3) / no_of_el;
                    end
                end
            end
            %             I_all_prev = I_all;
            %             I_all = zeros(no_of_radii,1,size(I_all,3));
            %             for (ind1=1:no_of_radii)
            %                 for (ind3=1:size(I_all,3))
            %                     I_now = I_all_prev(ind1,seg_range_use,ind3);
            %                     I_all(ind1,1,ind3) = mean(I_now(I_now >= 0));
            %                 end
            %             end
        else
            if (~isempty(seg_range))
                % remove the not needed segments
                no_of_segments = length(seg_range);
                I_all_prev = I_all;
                I_all = zeros(no_of_radii,no_of_segments,no_of_points);
                for (ind2=1:no_of_segments)
                    I_all(:,ind2,:) = I_all_prev(:,seg_range(ind2),:);
                end
            end
        end
        
        % average over radii, if specified
        if (plot_angle)
            no_of_radii = 1;
            if (isempty(radius_range))
                radius_range_use = 1:size(I_all,1);
            else
                radius_range_use = radius_range;
            end
            
            % If averaging is performed, i.e., more than one intensity
            % value is available, then the number of pixels in each
            % integration area needs to be known from the integration masks
            % data.
            if (isempty(integ_data))
                integ_data.integ_masks.norm_sum = ones(size(I_all,1),no_of_segments);
                if (length(radius_range_use) > 1)
                    fprintf('Warning: The integration mask could not be loaded from the file %s.\n',...
                        filename_integ_masks);
                    fprintf('Therefore the averaging over radii is not normalized by the number of pixels, which strongly influences the result.\n');
                end
            end
            
            % average over all pixels with positive intensities, i.e., skip
            % negative intensities
            I_all_prev = I_all;
            I_all = zeros(1,no_of_segments,no_of_points);
            for (ind2=1:no_of_segments)
                for (ind3=1:no_of_points)
                    no_of_pixels = 0;
                    for (ind1=1:length(radius_range_use))
                        if (I_all_prev(radius_range_use(ind1),ind2,ind3) >= 0)
                            I_all(1,ind2,ind3) = I_all(1,ind2,ind3) + ...
                                integ_data.integ_masks.norm_sum(radius_range_use(ind1),ind2) * I_all_prev(radius_range_use(ind1),ind2,ind3);
                            no_of_pixels = no_of_pixels + integ_data.integ_masks.norm_sum(radius_range_use(ind1),ind2);
                        end
                    end
                    if (no_of_pixels > 0)
                        I_all(1,ind2,ind3) = I_all(1,ind2,ind3) / no_of_pixels;
                    end
                end
            end
        end
        
        % average over points, if specified
        if (point_avg)
            if (isempty(point_range))
                point_range_use = 1:size(I_all,3);
            else
                point_range_use = point_range;
            end
            % average over all pixels with positive intensities, i.e., skip
            % negative intensities
            I_all_prev = I_all;
            I_all = zeros(no_of_radii,size(I_all,2),1);
            for (ind1=1:no_of_radii)
                for (ind2=1:no_of_segments)
                    no_of_el = 0;
                    for (ind3=1:length(point_range_use))
                        if (I_all_prev(ind1,ind2,point_range_use(ind3)) >= 0)
                            I_all(ind1,ind2,1) = I_all(ind1,ind2,1) + I_all_prev(ind1,ind2,point_range_use(ind3));
                            no_of_el = no_of_el +1;
                        end
                    end
                    if (no_of_el > 1)
                        I_all(ind1,ind2,1) = I_all(ind1,ind2,1) / no_of_el;
                    end
                end
            end
            
            no_of_points = 1;
        else
            if (~isempty(point_range))
                no_of_points = length(point_range);
            else
                no_of_points = size(I_all,3);
            end
        end
        
        if (~new_fig)
            legend_str = cell(1,no_of_segments);
        end
        
        % calculate q:
        if (isfield(d,'q')) && (~isempty(d.q))
            % define q_A aalways, if available
            q_A = d.q;
        else
            % calculate q, if it is not available in the (historic) data set,
            % and if the necessary parameters are provided
            if (plot_q || (~isempty(q_mul_pow)))
                error_base_str1 = 'You requested PlotQ or QMulPow, but the q-vector does not exist in the file. Please provide ';
                error_base_str2 = ' in order to calculate the q-vector.';
                if isempty(E_keV)
                    error('%s''E_keV''%s',error_base_str1,error_base_str2)
                end
                if isempty(det_dist_mm)
                    error('%s''DetDist_mm''%s',error_base_str1,error_base_str2)
                end
                if isempty(pixel_size_mm)
                    error('%s''PixelSize_mm''%s',error_base_str1,error_base_str2)
                end
                if isempty(radius)
                    error('You requested PlotQ or QMulPow, but neither the q-vector nor radius values are included in the file and thus q cannot be computed.');
                end
                q_A = pixel_to_q( radius, pixel_size_mm, det_dist_mm, E_keV );
            else
                % q_A is not needed -- set it to empty to indicate this
                q_A = [];
            end
        end
        % set some x-axis related values for plotting
        if (plot_angle)
            no_of_x_values = size(I_all,2);
            if (no_of_x_values < 2)
                error('At least two segments must be present for a plot as the function of the azimuthal angle.\n');
            end
            % calculate the azimuthal axis values for the angular plot
            if isfield(d,'phi_det')
                x_values = d.phi_det;
            else    % For backwards compatilibity with integrated data without angle
                x_values = 0:(360/no_of_x_values):(360-0.99*360/no_of_x_values);
            end
            % corresponding axis label
            x_label = '\Theta [ ^\circ ]';
        else
            if (plot_q)
                x_values = q_A;
                if (inverse_nm)
                    x_values = x_values * 10;
                    x_label = 'q [ nm^{-1} ]';
                else
                    x_label = 'q [ A^{-1} ]';
                end
            else
                if not((isfield(d,'radius')))
                    x_values = d.q;
                    x_label = 'q [ A^{-1} ]';
                else
                    x_values = radius;
                    x_label = 'pixel no.';
                    
                end
            end
        end
        
        % exclude zero or negative radii in case of logarithmic x-scale
        if (xlog)
            ind_x = find(x_values > 0);
        else
            ind_x = 1:length(x_values);
        end
        
        % set intensity multiplication values
        I_times = ones(length(x_values),1,1);
        y_label = 'average counts per pixel';
        if (~isempty(q_mul_pow))
            I_times(:,1,1) = q_A .^ q_mul_pow;
            y_label = [ y_label ' \times (q [ A^{-1} ])^{' ...
                num2str(q_mul_pow,'%.1f') '}' ];
        end
        
        if (isempty(y_values_returned))
            d1 = 0;
            d2 = 0;
            d3 = 0;
        else
            d1 = size(y_values_returned,1);
            d2 = size(y_values_returned,2);
            d3 = size(y_values_returned,3);
        end
        
        for (ind_point = 1:no_of_points)
            % determine the number of the point to be plotted
            if ((point_avg) || (isempty(point_range)))
                plot_ind_point = ind_point;
            else
                plot_ind_point = point_range(ind_point);
            end
            
            if (plot_angle)
                % plot as a function of the azimuthal angle
                
                % check for negative y-values
                if (ylog)
                    ind_y = find( I_all(1,:,plot_ind_point) > 0 );
                else
                    ind_y = 1:no_of_segments;
                end
                ind = intersect(ind_x,ind_y);
                
                % plot the intensity as a function of the azimuthal angle
                x_values_plotted = x_values(ind);
                y_values_plotted = I_all(1,ind,plot_ind_point);
                if (show_fig)
                if (fig_no > 0)
                    plot_function(x_values_plotted, y_values_plotted);
                end
                end
                
                % store the plotted values in the return array (untested)
                y_values_returned((d1+1):(d1+1), (d2+1):(d2+length(ind)), (d3+ind_point):(d3+ind_point)) = ...
                    y_values_plotted;
                 x_values_returned((d1+1):(d1+1), (d2+1):(d2+length(ind)), (d3+ind_point):(d3+ind_point)) = ...
                    x_values_plotted;               
                % do not add a legend
                legend_mul_seg = 0 ;
                hold all;
            else
                % plot as a function of radius or q
                
                for (ind_seg = 1:no_of_segments)
                    % determine the number of the segment to be plotted
                    if ((seg_avg) || (isempty(seg_range)))
                        plot_ind_seg = ind_seg;
                    else
                        plot_ind_seg = seg_range(ind_seg);
                    end
                    
                    % check for negative y-values
                    if (ylog)
                        ind_y = find( I_all(:, ind_seg,plot_ind_point) > 0 );
                    else
                        ind_y = 1:size(I_all,1);
                    end
                    ind = intersect(ind_x,ind_y);
                    
                    % plot the segment
                    x_values_plotted = x_values(ind);
                    y_values_plotted = I_all(ind, ind_seg,plot_ind_point) .* I_times(ind,1,1);
                    if show_fig
                    if (fig_no > 0)
                        plot_function(x_values_plotted, y_values_plotted);
                    end
                    end
                    % store the plotted values in the return array
                    y_values_returned((d1+1):(d1+length(ind)), (d2+ind_seg):(d2+ind_seg), (d3+ind_point):(d3+ind_point)) = ...
                        y_values_plotted;
                    x_values_returned((d1+1):(d1+length(ind)), (d2+ind_seg):(d2+ind_seg), (d3+ind_point):(d3+ind_point)) = ...
                        x_values_plotted;                    
                    % add the segment to the legend
                    if (fig_no > 0)
                        if (~new_fig)
                            legend_str{ind_seg} = [ 'seg. ' num2str(plot_ind_seg,'%03d') ];
                        else
                            legend_str(file_ind,:) = ...
                                [ num2str(file_ind,'%04d') ' seg.' num2str(plot_ind_seg,'%03d') ];
                        end
                        hold all;
                    end
                end
            end
            if (show_fig)
            if (fig_no > 0)
                if (clear_fig > 1)
                    hold off;
                else
                    hold all;
                end
            end
            title_str = strrep(filename,'\','\\');
            title_str = strrep(title_str,'_','\_');
            if (no_of_points > 1)
                title_str = [ title_str ', point ' num2str(plot_ind_point-1) ];
            end
            if (fig_no > 0)
                title( title_str );
                axis tight;
                xlabel(x_label);
                ylabel(y_label);
                if (~isempty(axis_scale))
                    axis( axis_scale );
                end
                if ((legend_mul_seg) && (no_of_segments > 1))
                    legend(char(legend_str));
                end
                
                drawnow;
            end
            first_plot = 0;
            
            if (sleep_time > 0.0)
                pause(sleep_time);
            end
        end
        
        if ((new_fig) && (fig_no > 0))
            fig_no = fig_no +1;
        end
        end
        
    end
    %     if ((~new_fig) && (size(legend_str,1) <= 10))
    %         legend(char(legend_str));
    %     end
end


function [I_all] = calc_I_point_avg(I_all_prev,point_range)

% average over all pixels with positive intensities, i.e., skip
% negative intensities
no_of_radii = size(I_all_prev,1);
no_of_segments = size(I_all_prev,2);
I_all = zeros(no_of_radii,no_of_segments,1);
for (ind1=1:no_of_radii)
    for (ind2=1:no_of_segments)
        no_of_el = 0;
        for (ind3=1:length(point_range))
            if (I_all_prev(ind1,ind2,point_range(ind3)) >= 0)
                I_all(ind1,ind2,1) = I_all(ind1,ind2,1) + I_all_prev(ind1,ind2,point_range(ind3));
                no_of_el = no_of_el +1;
            end
        end
        if (no_of_el > 1)
            I_all(ind1,ind2,1) = I_all(ind1,ind2,1) / no_of_el;
        end
    end
end
