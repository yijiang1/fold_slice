% Call function without arguments for instructions on how to use it

% Filename: $RCSfile: tune_valid_mask.m,v $
%
% $Revision: 1.5 $  $Date: 2012/09/02 15:13:40 $
% $Author: bunk $
% $Tag: $
%
% Description:
% remove outlyers of intensity that deviates from the azimuthal integration
% from the valid pixel mask
%
% Note:
% Call without arguments for a brief help text.
%
% Dependencies: 
% - image_read
%
% history:
%
% May 21st 2010, Oliver Bunk:
% 1st version

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

function [valid_mask] = tune_valid_mask(data_dir, varargin)
import beamline.radial_integ
import io.image_read
import plotting.display_valid_mask
import utils.find_files


% set default values for the variable input arguments:
% directory with the integrated data files
indir_integ_data = '~/Data10/analysis/radial_integration/';
% filename of the integrated data, empty to determine it from the raw data
% file name
filename_integ_data = [];
% use all cbf files
filename_mask = '*.cbf';
% filename for loading and saving the valid pixel mask
filename_valid_mask = '~/Data10/analysis/data/pilatus_valid_mask.mat';
% integration masks
filename_integ_masks = '~/Data10/analysis/data/pilatus_integration_masks.mat';

% size of the median filter that is use to smooth the data for identifying
% outlyers
median_size = 11;
% first pixel to start at
radius_from = 20;
% last pixel to check
radius_to = 0;
% only intensities above this threshold are considered for being hot
threshold_hot = 5;
% this value times the standard deviation of the intensity is used as hot pixel
% threshold
threshold_median = 3.0;
% save the updated mask
save_data = 0;
% display result in this figure
fig_no = 201;
% matching files to use
point_range = [];

% check minimum number of input arguments
if (nargin < 1)
    fprintf('\nUsage:\n');
    fprintf('[valid_mask]=%s(data_dir [[,<name>,<value>]...]);\n',mfilename);
    fprintf('Remove outlyers from the valid pixel mask by comparing azimuthally integrated data\n');
    fprintf('against the same data median filtered and rejecting pixels with a deviation\n');
    fprintf('in intensity specified in multiples of the standard deviation.\n');
    fprintf('\n');
    fprintf('The optional <name>,<value> pairs are:\n');
    fprintf('''FilenameMask'',<file specifier>      specify the files to be used from the data directory, empty string for all, default is ''%s''\n',...
        filename_mask);
    fprintf('''PointRange'',<vector or []>          matching files to use, default is [] for all files\n');
    fprintf('''IndirIntegData'',<filename.mat>      directory with the azimuthally integrated data, default is %s\n',...
        indir_integ_data);
    fprintf('''FilenameIntegData'',<filename.mat>   filename for the azimuthally integrated data, empty to determine the name\n');
    fprintf('                                     from the first raw data file name, default is ''%s''\n',...
        filename_integ_data);
    fprintf('''FilenameIntegMasks'',<filename>      Matlab file containing the integration masks, default is ''%s''\n',filename_integ_masks);
    fprintf('''RadiusFrom'',<integer>               no. of the pixel to start with, default is %.0f\n',radius_from);
    fprintf('''RadiusTo'',<integer>                 no. of the last pixel to check, default is %.0f\n',radius_to);
    fprintf('''MedianSize'',<integer>               size of the median filter in pixels, default is %.0f\n',...
        median_size);
    fprintf('''ThresholdHot'',<float>               pixels above this value are considered for being hot, default is %d\n',...
        threshold_hot);
    fprintf('''ThresholdMedian'',<float>            pixels outside the range (I+/-threshold_median*sqrt(I))\n');
    fprintf('                                     of the median filtered data are considered to be hot,\n');
    fprintf('                                     default is %.1f\n',...
        threshold_median);
    fprintf('''SaveData'',<0-no,1-yes>              save the valid pixel mask, default is %d\n',save_data);
    fprintf('''FilenameValidMask'',<path and filename>  Matlab file with the valid pixel indices,\n');
    fprintf('                                     default is %s\n',filename_valid_mask);
    fprintf('''FigNo'',<integer>                    number of the figure in which the result is displayed, default is %d\n',...
        fig_no);
    fprintf('\n');
    fprintf('Examples:\n');
    fprintf('[valid_mask]=%s(''~/Data10/pilatus/S05000-05999/S05715/e12612_1_05715_00000_00000.cbf'');\n',...
        mfilename);
    fprintf('[valid_mask]=%s(''~/Data10/pilatus/S05000-05999/S05715/*.cbf'');\n',...
        mfilename);

    error('At least the filename of the raw data has to be specified as input parameter.');
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
            no_of_in_arg = no_of_in_arg -1 + length(varargin);
        end
    end
end

% check number of input arguments
if (rem(no_of_in_arg,2) ~= 1)
    error('The optional parameters have to be specified as ''name'',''value'' pairs');
end


% parse the variable input arguments:
% initialize the list of unhandled parameters
vararg_remain = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'FilenameIntegMasks' 
            filename_integ_masks = value;
        case 'FilenameMask'
            filename_mask = value;
        case 'PointRange'
            point_range = value;
        case 'IndirIntegData'
            indir_integ_data = value;
        case 'FilenameIntegData'
            filename_integ_data = value;
        case 'RadiusFrom',
            radius_from = round(value);
        case 'RadiusTo',
            radius_to = round(value);
        case 'MedianSize'
            median_size = round(value);
        case 'ThresholdMedian'
            threshold_median = value;
        case 'ThresholdHot'
            threshold_hot = value;
        case 'FilenameValidMask' 
            filename_valid_mask = value;
        case 'SaveData' 
            save_data = value;
        case 'FigNo' 
            fig_no = value;
        otherwise
            vararg_remain{end+1} = name; %#ok<AGROW>
            vararg_remain{end+1} = value; %#ok<AGROW>
    end
end

vararg_remain{end+1} = 'UnhandledParError';
vararg_remain{end+1} = 0;
vararg_remain{end+1} = 'DisplayFilename';
vararg_remain{end+1} = 0;

% set some default values for the plot window
set(0, 'DefaultAxesfontsize', 12);
set(0, 'DefaultAxeslinewidth', 1, 'DefaultAxesfontsize', 12);
set(0, 'DefaultLinelinewidth', 1);

% get all matching filenames
if (data_dir(end) ~= '/')
    data_dir(end+1) = '/';
end
[data_dir,fnames,vararg_remain] = ...
    find_files( [ data_dir filename_mask ], vararg_remain );

if (length(fnames) < 1)
    error('No matching files found for %s%s.\n',data_dir,filename_mask);
end

% load the current valid pixel mask in variable valid_mask
fprintf('loading the existing valid mask %s\n',filename_valid_mask);
load(filename_valid_mask);
framesize = valid_mask.framesize(1) * valid_mask.framesize(2);

% load the integration masks in variable integ_masks
fprintf('Loading the integration masks from %s\n',filename_integ_masks);
load(filename_integ_masks);
no_of_radii = length(integ_masks.radius);

if ((radius_to < radius_from) || (radius_to > no_of_radii))
    radius_to = no_of_radii;
end

% process the frames
ind_hot = [];
ind_dark = [];
integ_data = [];
fprintf('data directory is %s\n',data_dir);
if (isempty(point_range))
    point_range = 1:length(fnames);
else
    ind = find(point_range <= length(fnames));
    if (length(point_range) ~= length(ind))
        fprintf('Warning, %d value(s) from the specified point range are out of the range [1,%.0f] and not used.\n',...
            length(point_range)-length(ind),length(fnames));
        point_range = point_range(ind);
    end
end

for (point_ind=1:length(point_range)) 
    f_ind = point_range(point_ind);
    % read the raw data
    fprintf('%3d/%3d: reading %s%s\n',f_ind,length(point_range),...
        data_dir,fnames(f_ind).name);
    filename_raw = [data_dir fnames(f_ind).name ];
    [frame] = image_read(filename_raw,vararg_remain);

    % check that the files have identical dimensions
    if ((size(frame.data,1) ~= valid_mask.framesize(1)) || ...
        (size(frame.data,2) ~= valid_mask.framesize(2)))
            error('The valid pixel mask has %d x %d pixels, this frame has %d x %d pixels',...
            valid_mask.framesize(1),valid_mask.framesize(2),...
            size(frame.data,1),size(frame.data,2));
    end

    % read the radially integrated data
    if (isempty(integ_data))
        % determine filename for the integrated data from the first raw
        % data filename
        if (isempty(filename_integ_data))
            [pathstr, filename_integ_data] = fileparts(fnames(f_ind).name);
            filename_integ_data = [ filename_integ_data '_integ.mat' ]; %#ok<AGROW>
        end
        filename_integ_data = fullfile(indir_integ_data,filename_integ_data);
        fprintf('Loading the integrated intensities from %s\n',...
            filename_integ_data);
        integ_data = load(filename_integ_data);
        
        % take the median of all segments with positive intensities, i.e.,
        % skip negative intensities
        I_all_prev = integ_data.I_all;
        no_of_segments = size(I_all_prev,2);
        no_of_points = size(I_all_prev,3);
        I_all = zeros(no_of_radii,no_of_points);
        I_std = zeros(no_of_radii,no_of_points);
        if (no_of_segments > 1)
            fprintf('Using the median of %d segments.\n',no_of_segments);
        end
        for (ind1=1:no_of_radii)
            for (ind3=1:no_of_points)
                no_of_el = 0;
                I_use = zeros(1,no_of_segments);
                ind_I_use = zeros(1,no_of_segments);
                for (ind2=1:no_of_segments)
                    if (I_all_prev(ind1,ind2,ind3) >= 0)
                        no_of_el = no_of_el +1;
                        I_use(no_of_el) = I_all_prev(ind1,ind2,ind3);
                        ind_I_use(no_of_el) = ind2;
                    end
                end
                if (no_of_el > 1)
                    [I_sorted,ind_sorted] = sort(I_use(1:no_of_el));
                    ind_median = round(0.5*no_of_el);
                    I_all(ind1,ind3) = I_sorted(ind_median);
                    % get the standard deviation of this intensity
                    I_std(ind1,ind3) = integ_data.I_std(ind1,ind_I_use(ind_sorted(ind_median)),ind3);
                end
            end
        end        
        
        % print this information once rather than for each file
        fprintf('Checking radii from %d to %d.\n',radius_from,radius_to);
    end

    % get the index to the integrated data
    ind = 1;
    ind_max = length(integ_data.filenames_all);
    while ((ind <= ind_max) && ...
            (isempty(strfind(integ_data.filenames_all{ind},fnames(f_ind).name))))
        ind = ind +1;
    end
    if (ind > ind_max)
        error('Could not find integrated data for raw data file %s in %s.',...
            filename_raw,filename_integ_data);
    end
    data_integ = squeeze(I_all(:,ind));
    data_integ_std = squeeze(I_std(:,ind));
    
    % median filtered data for comparison
    data_integ_med = medfilt1(data_integ,median_size,size(data_integ,1),1);

%     figure(fig_no+2);
%     hold off;
%     clf;
%     semilogy(data_integ);
%     hold all;
%     semilogy(data_integ_med);
%     semilogy(data_integ_med+data_integ_std*threshold_median);
%     semilogy(data_integ_med-data_integ_std*threshold_median);

    frame_cmp = ones(valid_mask.framesize) -2;
    frame_cmp_std = zeros(valid_mask.framesize);
    for (ind_r = radius_from:radius_to)
        for (ind_seg = 1:no_of_segments)
            if (integ_masks.norm_sum(ind_r,ind_seg) > 0)
                frame_cmp(integ_masks.indices{ind_r,ind_seg}) = ...
                    data_integ_med(ind_r);
                frame_cmp_std(integ_masks.indices{ind_r,ind_seg}) = ...
                    data_integ_std(ind_r);
            end
        end
    end

    % 
    ind_dark = union(ind_dark, ...
        find((frame.data >= 0) & ...
            (frame_cmp >= 0) & ...
            (frame.data < frame_cmp - threshold_median*frame_cmp_std)));
    % only consider pixels of sufficient intensity for being hot
    ind_hot = union(ind_hot, ...
        find((frame.data > threshold_hot) & ...
            (frame_cmp >= 0) & ...
            (frame.data > frame_cmp + threshold_median*frame_cmp_std)));
end


% calculate the complementary masks of the valid pixels
valid_mask.indices = intersect(valid_mask.indices,...
    setdiff(1:framesize,union(ind_dark,ind_hot)));

fprintf('In total %d dark and %d hot pixels found.\n',...
    length(ind_dark),length(ind_hot));
fprintf('%d valid pixels remain.\n',length(valid_mask.indices));


if (save_data)
    % create a backup of the mask
    if (exist(filename_valid_mask,'file'))
        filename_mask_backup = [ filename_valid_mask '.bak' ];
        fprintf('Copying the current mask %s to %s\n',filename_valid_mask,...
            filename_mask_backup);
        copyfile(filename_valid_mask,filename_mask_backup);
    end

    % save the masks
    fprintf('Saving valid_mask to %s\n',filename_valid_mask);
    save(filename_valid_mask,'valid_mask');
    % plot new valid pixel mask
    display_valid_mask('FilenameValidMask',filename_valid_mask,...
        'NoHelp',1,'FigNo',fig_no);
else
    fprintf('The updated valid pixel mask is NOT saved.\n');
end

% plot the additional invalid pixels
figure(fig_no+1);

% mark the valid pixels as 1, leave the invalid at 0
frame = zeros(valid_mask.framesize);
frame(valid_mask.indices) = 1;
frame(ind_dark) = -10;
frame(ind_hot) = 10;
imagesc(frame);
caxis([-10 10]);
axis xy;
axis equal;
axis tight;
colorbar;
title_str = ['valid pixels, ' ...
    num2str(length(ind_dark)+length(ind_hot),'%d') ...
    ' update(s) marked with intensity -10/10'];
title(title_str);
set(gcf,'Name','valid pixels, updates marked');
