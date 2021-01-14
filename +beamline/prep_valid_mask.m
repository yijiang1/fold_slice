% Call function without arguments for instructions on how to use it

% Filename: $RCSfile: prep_valid_mask.m,v $
%
% $Revision: 1.8 $  $Date: 2016/01/21 15:07:41 $
% $Author: guizar_m $
% $Tag: $
%
% Description:
% prepare a list of the linear indices for the valid pixels
%
% Note:
% Call without arguments for a brief help text.
%
% Dependencies: 
% - image_read
%
% history:
%
% May 15th 2010, Oliver Bunk:
% add command line argument for ThresholdMedian
%
% September 4th 2009, Oliver Bunk: 
% use find_files rather than dir to find the files
%
% May 9th 2008, Oliver Bunk: 1st documented version

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


function [valid_mask] = prep_valid_mask(data_dir, varargin)
import io.image_read
import plotting.display_valid_mask
import utils.find_files

% initialize return arguments
valid_mask = struct('indices',[], 'framesize',[]);

% set default values for the variable input arguments:
% use all cbf files
filename_mask = '*.cbf';
% filename for loading and saving the valid pixel mask
filename_valid_mask = '~/Data10/analysis/data/pilatus_valid_mask.mat';
% below this threshold intensity a pixel is considered to be dark
threshold_dark = 1;
% above this threshold intensity a pixel is considered to be hot
threshold_hot = 20;
% this value times the square root of the intensity is used as hot pixel
% threshold
threshold_median = 5.0;
% replace the existing mask
extend = 'no';
% save the mask
save_data = 1;
% display result in this figure
fig_no = 200;

% check minimum number of input arguments
if (nargin < 1)
    fprintf('\nUsage:\n');
    fprintf('[valid_mask]=%s(data_dir [[,<name>,<value>]...]);\n',mfilename);
    fprintf('Prepare a list of the linear indices for the valid pixels.\n');
    fprintf('To get reliable data a series of at least 10 frames should be analyzed.\n');
    fprintf('The direct beam region will be regarded as invalid since it is out of the\n');
    fprintf('range for valid pixels. To ''repair'' this one should take a second series of\n');
    fprintf('exposures at a different detector position and call this macro with the ''Extend'',''or''\n');
    fprintf('option.\n');
    fprintf('\n');
    fprintf('The optional <name>,<value> pairs are:\n');
    fprintf('''FilenameMask'',<file specifier>     specify the files to be used from the data directory, empty string for all, default is ''%s''\n',...
        filename_mask);
    fprintf('''ThresholdDark'',<float>             pixels permanently below this value are considered to be dark, default is %d\n',...
        threshold_dark);
    fprintf('''ThresholdHot'',<float>              pixels at least once above this value are considered to be hot, default is %d\n',...
        threshold_hot);
    fprintf('''ThresholdMedian'',<float>           pixels of intensity I above the constant ThresholdHot and above\n');
    fprintf('                                    ThresholdMedian times (I+sqrt(I)) are considered to be hot, 0 to deactivate this additional threshold,\n');
    fprintf('                                    default is %.1f\n',...
        threshold_median);
    fprintf('''SaveData'',<0-no,1-yes>             save the valid pixel mask, default is %d\n',save_data);
    fprintf('''FilenameValidMask'',<path and filename>  Matlab file with the valid pixel indices,\n');
    fprintf('                                    default is %s\n',filename_valid_mask);
    fprintf('''Extend'',<''and'', ''or'' or ''no''>  update an existing mask using the specified conjunction, default is %s\n',...
        extend);
    fprintf('''FigNo'',<integer>                   number of the figure in which the result is displayed, default is %d\n',...
        fig_no);
    fprintf('\n');
    fprintf('Examples:\n');
    fprintf('[valid_mask]=%s(''~/Data10/pilatus/air_scattering/'');\n',...
        mfilename);
    fprintf('[valid_mask]=%s(''~/Data10/pilatus/air_scattering_det_pos_2/'',''Extend'',''or'');\n',...
        mfilename);

    error('At least the data directory has to be specified as input parameter.');
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
        case 'ThresholdDark'
            threshold_dark = value;
        case 'ThresholdHot'
            threshold_hot = value;
        case 'ThresholdMedian'
            threshold_median = value;
        case 'FilenameMask' 
            filename_mask = value;
        case 'FilenameValidMask' 
            filename_valid_mask = value;
        case 'SaveData' 
            save_data = value;
        case 'FigNo' 
            fig_no = value;
        case 'Extend' 
            extend = value;
        otherwise
            vararg_remain{end+1} = name; %#ok<AGROW>
            vararg_remain{end+1} = value; %#ok<AGROW>
    end
end

vararg_remain{end+1} = 'UnhandledParError';
vararg_remain{end+1} = 0;

% check extend parameter
if ((~strcmp(extend,'no')) && ...
    (~strcmp(extend,'and')) && (~strcmp(extend,'or')))
    error('extend must be ''and'', ''or'' or ''no''\n');
end

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

if (~strcmp(extend,'no'))
  if exist(filename_valid_mask,'file')
    fprintf('loading the existing valid mask %s\n', ...
      filename_valid_mask);
    load(filename_valid_mask);
    ind_existing_valid = valid_mask.indices;
  else
    fprintf('no prior valid mask %s found\n', ...
      filename_valid_mask);
    ind_existing_valid = '';
  end
end


% process the frames
ind_hot = [];
ind_dark = [];
fprintf('data directory is %s\n',data_dir);
for (f_ind=1:length(fnames)) 
    fprintf('%3d/%3d: reading %s%s\n',f_ind,length(fnames),...
        data_dir,fnames(f_ind).name);
    [frame] = image_read([data_dir fnames(f_ind).name ],vararg_remain);
    frame.data = double(frame.data);
    for (frame_ind = 1:size(frame.data,3))
        % median filtered data for comparison
        if (threshold_median ~= 0)
            data_med = frame.data(:,:,frame_ind);

            % add pixels at the module boundary to ease median filtering
            ind = find(data_med == 0);
            data_med_shift = circshift(data_med,[2 2]);
            data_med(ind) = data_med_shift(ind);

            ind = find(data_med == 0);
            data_med_shift = circshift(data_med,[-2 -2]);
            data_med(ind) = data_med_shift(ind);

            ind = find(data_med == 0);
            data_med_shift = circshift(data_med,[-2 2]);
            data_med(ind) = data_med_shift(ind);

            ind = find(data_med == 0);
            data_med_shift = circshift(data_med,[2 -2]);
            data_med(ind) = data_med_shift(ind);

            % median filter the data
            data_med = medfilt2(data_med,[5 5]);

            % the square root of the intensity estimates the standard deviation
            data_med_sqrt = data_med.^0.5;
        end

        if (f_ind == 1)
            framesize1 = size(frame.data,1);
            framesize2 = size(frame.data,2);
            framesize = framesize1 * framesize2;
        end

        % check that the file have identical dimensions
        if ((framesize1 ~= size(frame.data,1)) || ...
            (framesize2 ~= size(frame.data,2)))
                error('The previous file(s) had %d x %d pixels, this frame has %d x %d pixels',...
                framesize1,framesize2,size(frame.data,1),size(frame.data,2));
        end

        % pixels are considered to be dark if the intensity is below the
        % constant threshold
        ind = find(frame.data(:,:,frame_ind) < threshold_dark);
        fprintf('%6d dark pixels below %10.3e counts, ', ...
            length(ind),threshold_dark);
        if (f_ind == 1)
            ind_dark = ind;
        else
            % dark pixels must be dark in all frames
            ind_dark = intersect(ind_dark,ind);
        end

        % hot pixels are hot if they are above the threshold
        ind = find(frame.data(:,:,frame_ind) > threshold_hot);
        % and, if active, above the intensity plus a threshold times the square
        % root of the intensity as an estimation of the countin statistics
        % error
        if (threshold_median ~= 0.0)
            ind = intersect(ind,find((frame.data(:,:,frame_ind) > data_med+threshold_median*data_med_sqrt)));
            fprintf('%4d hot pixels above %d and %.1f * sqrt(intensity) counts\n', ...
                length(ind),threshold_hot,threshold_median);
        else
            fprintf('%4d hot pixels above %d counts\n', ...
                length(ind),threshold_hot);
        end
        % for hot pixels it is enough to be above the threshold in one frame
        ind_hot = union(ind_hot,ind);   
    end
end

% calculate the complementary masks of the valid pixels
valid_mask.indices = setdiff(1:framesize,union(ind_dark,ind_hot));

fprintf('In total %d dark and %d hot pixels found.\n',...
    length(ind_dark),length(ind_hot));
fprintf('%d valid pixels remain.\n',length(valid_mask.indices));

if (~strcmp(extend,'no'))
     fprintf('Extending the existing valid pixel mask of %d pixels\n',...
             length(ind_existing_valid));
     if (strcmp(extend,'and'))
         fprintf('using the and conjugation\n');
         valid_mask.indices = ...
             intersect(valid_mask.indices,ind_existing_valid);
     else
         fprintf('using the or conjugation\n');
         if ~isempty(ind_existing_valid)
             valid_mask.indices = ...
                 union(valid_mask.indices,ind_existing_valid);
         end
     end
     fprintf('The combined mask has %d valid pixels.\n',...
             length(valid_mask.indices));
end

% store the frame size in the return data
valid_mask.framesize = [framesize1 framesize2];

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
end

% plot new valid pixel mask
if (fig_no > 0)
  display_valid_mask('FilenameValidMask',filename_valid_mask,...
    'NoHelp',1,'FigNo',fig_no);
end
