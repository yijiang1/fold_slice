% Call function without arguments for instructions on how to use it

% Filename: $RCSfile: radial_integ.m,v $
%
% $Revision: 1.12 $  $Date: 2016/01/21 15:11:50 $
% $Author:  $
% $Tag: $
%
% Description:
% radial integration of 2D data read from file(s)
%
% Note:
% Call without arguments for a brief help text.
% The integration masks need to be prepared first using prep_integ_masks.m
%
% Dependencies:
% - image_read
%
% history:
%
% February 18 2015:
% updated to use new function names of parallel toolbox in Matlab 2014b
%
% July 22nd 2010:
% add simple parallel processing using parfor
%
% April 28th 2010:
% use default_parameter_value
%
% June 5th 2008: 1st documented version

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

function [ I,vararg_remain ] = radial_integ(filename_masks,varargin)
import beamline.prep_integ_masks
import io.image_read
import plotting.plot_radial_integ
import utils.default_parameter_value
import utils.find_files
import utils.abspath

% initialize return arguments
I = struct('I_all',[], 'I_std',[],'filenames_all',[],'q',[],'radius',[]);

% set default values for the variable input arguments:
outdir_data = default_parameter_value(mfilename,'OutdirData');
filename_integ_masks = default_parameter_value(mfilename,'FilenameIntegMasks');
r_max_forced = default_parameter_value(mfilename,'rMaxForced');
fig_no = default_parameter_value(mfilename,'FigNo');
save_combined_I = default_parameter_value(mfilename,'SaveCombinedI');
recursive = default_parameter_value(mfilename,'Recursive');
use_find = default_parameter_value(mfilename,'UseFind');
unhandled_par_error = default_parameter_value(mfilename,'UnhandledParError');
parallel_tasks_max = 1; %default_parameter_value(mfilename,'ParTasksMax');
save_format = '-v6';
use_mex = true;
c_reader = true;
useStack = true;

% check minimum number of input arguments
if (nargin < 1)
    fprintf('\nUsage:\n');
    fprintf('%s(filename_mask,  [[,<name>,<value>] ...]);\n',mfilename);
    fprintf('filename_mask can be something like ''*.cbf'' or ''image.cbf'' or\n');
    fprintf('a cell array of filenames or filename masks like {''dir1/*.cbf'',''dir2/*.cbf''}.\n');
    fprintf('The optional <name>,<value> pairs are:\n');
    fprintf('''OutdirData'',<directory>             save the integrated intensities to files in this directory, '''' for no saving, default is %s\n',outdir_data);
    fprintf('''FilenameIntegMasks'',<filename>      Matlab file containing the integration masks, default is ''%s''\n',filename_integ_masks);
    fprintf('''rMaxForced'',<radius in pixel>       stop integration at this maximum r even if the integration masks reach further, default is 0 - do not stop\n');
    fprintf('''FigNo'',<figure number>              number of the figure for an online plot of the intensities in case parallel processing is not used, 0 for no plot, default is %d\n',fig_no);
    fprintf('''SaveFormat'',<format string>   default is %s\n',save_format);
    fprintf('''SaveCombinedI'',<0-no, 1-yes>        save intensities from all specified files found in one directory in a single file, default is yes\n');
    fprintf('''Recursive'',<0-no, 1-yes>            recursively integrate files in all matching sub-directories, default is yes\n');
    fprintf('''ParTasksMax'',<integer>              specify the maximum number of CPU cores to use, 1 to deactivate the use of parallel computing, default is %d\n',parallel_tasks_max);
    fprintf('''UseFind'',<0-no, 1-yes>              use Linux/Unix command find to interprete the filename mask, default is yes\n');
    fprintf('''UseMex'', <0-no, 1-yes>              use radial_integ_mex; usually faster than MATLAB, default is yes\n');
    fprintf('''CReader'', <0-no, 1-yes>             use the fast measurement reader; usually  faster than image_read, default is yes\n');
    fprintf('''UseStack'', <0-no, 1-yes>            load all detector frames into memory before calling the radial_integ functions; default is yes\n');
    fprintf('''UnhandledParError'',<0-no,1-yes>     exit in case not all named parameters are used/known, default is %d\n',unhandled_par_error);
    fprintf('Examples:\n');
    fprintf('%s(''~/Data10/pilatus/mydatadir/*.cbf'',''OutdirData'',''~/Data10/analysis/radial_integ/'');\n',mfilename);
    fprintf('Additional <name>,<value> pairs recognized by image_read can be specified.\n');
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
vararg = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'OutdirData'
            outdir_data = value;
        case 'SaveFormat'
            save_format = value;
        case 'FilenameIntegMasks'
            filename_integ_masks = value;
        case 'rMaxForced'
            r_max_forced = value;
        case 'FigNo'
            fig_no = value;
        case 'SaveCombinedI'
            save_combined_I = value;
        case 'Recursive'
            recursive = value;
        case 'UseFind'
            use_find = value;
        case 'UnhandledParError'
            unhandled_par_error = value;
        case 'ParTasksMax'
            parallel_tasks_max = value;
        case 'UseMex'
            use_mex = value;
        case 'CReader'
            c_reader = value;
        case 'UseStack'
            useStack = value;
        otherwise
            vararg{end+1} = name; %#ok<AGROW>
            vararg{end+1} = value; %#ok<AGROW>
    end
end

% initialize the list of unhandled parameters
vararg_remain = cell(0,0);

% do not exit in image_par in case of unhandled parameters
if (~unhandled_par_error)
    vararg{end+1} = 'UnhandledParError';
    vararg{end+1} = 0;
end

if (~isempty(outdir_data))
    % add slash to output directory
    if (outdir_data(end) ~= '/')
        outdir_data = [ outdir_data '/' ];
    end
    
    % create output directory
    [mkdir_stat,mkdir_message] = mkdir(outdir_data);
    if (~mkdir_stat)
        error('invalid directory %s: %s',outdir_data,mkdir_message);
    end
    if ((mkdir_stat) && (isempty(mkdir_message)))
        fprintf('The output directory %s has been created.\n',outdir_data);
    else
        fprintf('The output directory is %s.\n',outdir_data);
    end
else
    fprintf('data are not saved\n');
end

% load integration masks from this file
% this loads:
% center_xy, no_of_segments, integ_masks
fprintf('loading integration masks from %s\n',filename_integ_masks);
load(filename_integ_masks);
if ((~exist('center_xy','var')) && (exist('center_x','var')))
    center_xy(1) = center_x;
    center_xy(2) = center_y;
    if (~exist('integ_masks','var'))
        integ_masks.radius = r;
        integ_masks.indices = masks_r;
        integ_masks.norm_sum = mask_r_sum;
    end
end
fprintf('center at (x, y) = (%.1f, %.1f)\n',center_xy(1),center_xy(2));

% limit radial range
if (r_max_forced > 0)
    ind = find( integ_masks.radius < r_max_forced );
    if (length(ind) < 1)
        fprintf('No radii below rMaxForced = %d found\n',r_max_forced);
        return;
    end
    integ_masks.radius = integ_masks.radius(1:ind(end));
    integ_masks.norm_sum = integ_masks.norm_sum(1:ind(end), :);
end
if isfield(integ_masks,'radius')
fprintf('radii from %d to %d\n',...
    integ_masks.radius(1),integ_masks.radius(end));
else
   fprintf('radii from %d to %d\n',...
    integ_masks.q(1),integ_masks.q(end));
end
% ease handling by ensuring that filename_masks is a cell array
if (~iscell(filename_masks))
    filename_masks = { filename_masks };
end


% initialize parallel processing if this is enabled and not yet done
if (parallel_tasks_max > 1)
    %matlabpool_size = matlabpool('size');
    %if (matlabpool_size < 1)
    if isempty(gcp('nocreate'))  %MGS2015  If there is no current pool
        % create a scheduler object using the default configuration, which is a
        % local scheduler if nothing else has been installed
        % scheduler = findResource('scheduler','type', defaultParallelConfig);
        scheduler = parcluster; %MGS2015
        
        % adapt maximum number of tasks/workers, if necessary
        %cluster_size = get(scheduler,'ClusterSize');
        cluster_size = scheduler.NumWorkers; %MGS2015
        if (parallel_tasks_max > cluster_size)
            fprintf('Adapting the maximum number of tasks from %d to %d.\n',...
                parallel_tasks_max, cluster_size);
            parallel_tasks_max = cluster_size;
        end
        
        % open a Matlab pool for simple parallel processing
        if (parallel_tasks_max > 1)
            %matlabpool('open',parallel_tasks_max);%MGS2015
            pool = parpool(parallel_tasks_max);
            fprintf('Using parallel processing with %d tasks.\n', ...
                parallel_tasks_max);
        end
    else
        pool = gcp;%MGS2015
        if ( pool.NumWorkers < parallel_tasks_max )
            fprintf('%s: usage of up to %d CPUs in parallel has been specified but an already open matlabpool with %d workers has been found and will be used instead\n', ...
                mfilename, parallel_tasks_max, pool.NumWorkers);
            parallel_tasks_max = pool.NumWorkers;
        end
    end
end

pool.IdleTimeout = Inf;

if ((parallel_tasks_max > 1) && (fig_no > 0))
    fprintf('%s: Online plotting is disabled since parallel processing is enabled.\n', ...
        mfilename);
end

% loop over all filename masks
ind_mask_max = length(filename_masks);

% Initialize variables for saving
no_of_segments = size(integ_masks.indices,2);
if isfield(integ_masks,'radius')
    radius = integ_masks.radius;
    ind_r_max = length(radius);
else
    radius = [];
    q = integ_masks.q;
    ind_r_max = length(q);
end
if isfield(integ_masks,'q')
    q = integ_masks.q;
else
    q = [];
end

for (ind_mask = 1:ind_mask_max) %#ok<*NO4LP>
    filename_mask = filename_masks{ind_mask};
    fprintf('%s:\n',filename_mask);
    [data_dir,fnames] = find_files( filename_mask, 'UseFind',use_find );
    
    if (length(fnames) < 1)
        fprintf('No matching files found for %s.\n',filename_mask);
        continue;
    end
    
    % collect recursively all matching file names
    [ filenames_all ] = ...
        collect_radial_integ_filenames(data_dir, fnames, ...
        recursive, ...
        vararg);
    
    % prepare for integration of the so far identified files
    file_ind_max = length(filenames_all);
    for file_ind=1:file_ind_max
        filenames_all{file_ind}=abspath(filenames_all{file_ind});
    end
    
    % get the number of frames per file by loading the first file (not very
    % elegant)
    [frame] = image_read(filenames_all{1}, vararg);
    no_of_frames = size(frame.data,3);
    I_all = zeros(ind_r_max, no_of_segments, no_of_frames, file_ind_max);
    I_std = zeros(ind_r_max, no_of_segments, no_of_frames, file_ind_max);

    if (parallel_tasks_max > 1)
        % integration using parallel processing
        parfor (file_ind = 1:file_ind_max)
            
            % read the raw data frame and integrate it
            [frame_I, frame_std] = ...
                perform_radial_integ_parallel(file_ind, file_ind_max, ...
                filenames_all{file_ind}, ...
                integ_masks, ind_r_max, no_of_segments, ...
                vararg);
            
            %             no_of_frames = size(frame_I,3);
            %             if (no_of_frames ~= size(I_all,3))
            %                 error('number of frames per file changes from %d to %d',size(I_all,3),no_of_frames);
            %             end
            
            I_all(:,:,:,file_ind) = frame_I;
            I_std(:,:,:,file_ind) = frame_std;
        end
    else
        % read the raw data
        if c_reader
            try
                [~, ~, ext] = fileparts(filenames_all{1});
                arg.data_path = filenames_all';
                arg.nthreads = min(round(feature('numcores')*0.8),14);
                arg.precision = 'single';
                arg.extension = ext(2:end);
                if strcmpi(ext(1:end), 'h5') && ~isempty(find(strcmp(varargin, 'H5Location')))
                    arg.data_location = varargin{find(strcmp(varargin, 'H5Location'))+1};
                end
                frameStorage.data = io.read_measurement(arg);
                frameStorage.data = permute(frameStorage.data,[2 1 3]);
                frameStorage.data = flip(flip(frameStorage.data,1),2);
            catch ME
                fprintf('Failed to load data. If the problem persists, set c_reader=false.\n');
                rethrow(ME);
            end
        else
            [frameStorage] = image_read(filename_masks, vararg);
        end
        if ~useStack
            for (file_ind = 1:file_ind_max)
                
                % read the raw data frame and integrate it
                [frame_I, frame_std] = ...
                    perform_radial_integ(file_ind, file_ind_max, ...
                    frameStorage.data(:,:,file_ind), ...
                    integ_masks, ind_r_max, no_of_segments, use_mex, ...
                    vararg);
                
                I_all(:,:,:,file_ind) = frame_I;
                I_std(:,:,:,file_ind) = frame_std;
                
                % plot integrated intensities as feedback
                if (fig_no > 0)
                    if isfield(integ_masks,'radius')
                        d.radius = radius;
                    else
                        d.radius= q;
                    end
                    d.I_all = frame_I;
                    d.I_std = frame_std;
                    plot_radial_integ(d,'FigNo',fig_no);
                    drawnow;
                end
            end
        else
                 % integrate it
                try
                    [frame_I, frame_std] = radial_integ_mex(int32(ind_r_max),int32(no_of_segments), integ_masks.norm_sum, integ_masks.indices, double(frameStorage.data));
                catch
                    tmpPath = fileparts(mfilename('fullpath'));
                    fprintf('Recompiling mex function...\n');
                    
                    % Fall back to single thread if the OpenMP fail.
                    
                    eval(['mex ' fullfile(tmpPath, 'private', 'radial_integ_mex.cpp') ' -outdir ' fullfile(tmpPath, 'private')]);
                    try
                        [frame_I, frame_std] = radial_integ_mex(int32(ind_r_max),int32(no_of_segments), integ_masks.norm_sum, integ_masks.indices, double(frameStorage.data));
                    catch ME
                        fprintf('radial_integ_mex failed. If the problem persists, consider setting use_mex=false.\n');
                        rethrow(ME);
                    end
                end
                I_all(:,:,1,:) = frame_I;
                I_std(:,:,1,:) = frame_std;           
        end
    end
    
    % reshuffle the data to get rid off the frame-within-file dimension,
    % dimension 3.
    % This would be easier with linear indexing in case the frame and file
    % dimensions would be 1 and 2.
    I_all_org = I_all;
    I_std_org = I_std;
    I_all = zeros(size(I_all_org,1), size(I_all_org,2), size(I_all_org,3) * size(I_all_org,4));
    I_std = zeros(size(I_all));
    for (ind_frame = 1:size(I_all_org,3))
        for (ind_file = 1:size(I_all_org,4))
            I_all(:,:,(ind_file-1)*size(I_all_org,3)+ind_frame) = I_all_org(:,:,ind_frame,ind_file);
            I_std(:,:,(ind_file-1)*size(I_std_org,3)+ind_frame) = I_std_org(:,:,ind_frame,ind_file);
        end
    end
    
    % save data, if this option is enabled
    if (~isempty(outdir_data))
        if (save_combined_I)
            % save all integrated frames as single Matlab file
            
            if  (exist('I_all','var'))
                % use first file as file-name base
                [~, name] = fileparts(filenames_all{1});
                %                name = name(1:end-12);
                fname_out = fullfile(outdir_data, [ name '_integ.mat' ]);
                fprintf('saving %s\n',fname_out);
                % remove directory information before storing the filenames
                for (file_ind = 1:file_ind_max)
                    [~, name, extension] = fileparts(filenames_all{file_ind});
                    filenames_all{file_ind} = [ name extension ];
                end
                norm_sum = integ_masks.norm_sum;
                save(fname_out,'I_all','I_std', 'norm_sum', 'filenames_all','radius','q','angular_segments','phi_det', save_format);
            else
                fprintf('No data to save for directory %s\n',data_dir);
            end
        else
            
            % save the integrated data for each frame as separate ASCII
            % file
            savedat = zeros(ind_r_max, no_of_segments +1);
            if isfield(integ_masks,'radius')
            savedat(:,1) = radius;
            else
                savedat(:,1)= q;
            end
            for (file_ind = 1:file_ind_max)
                % save integrated data for this image in the output arrays
                savedat(:,2:end) = I_all(:,:,file_ind);
                
                [pathstr, name] = fileparts(filenames_all{file_ind});
                fname_out = fullfile(pathstr, [ name '_integ.txt' ]);
                fprintf('saving %s\n',fname_out);
                save([outdir_data fname_out],'savedat','-ascii');
            end
            
            fprintf('\nOutput data format:\n');
            fprintf('- first column with radius of circle in pixel\n');
            fprintf('- further columns with average intensity in circle segment\n');
        end
    end
    
    % compile return value
    I(ind_mask).I_all = I_all;
    I(ind_mask).I_std = I_std;
    if isfield(integ_masks,'radius')
    I(ind_mask).radius = integ_masks.radius;
    end
    I(ind_mask).norm_sum = integ_masks.norm_sum;
    I(ind_mask).filenames_all = filenames_all;
    if isfield(integ_masks,'q')
    I(ind_mask).q = integ_masks.q;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ filenames_all ] = ...
    collect_radial_integ_filenames(data_dir, fnames, ...
    recursive, ...
    vararg)
import beamline.prep_integ_masks
import io.image_read
import plotting.plot_radial_integ
import utils.default_parameter_value
import utils.find_files
% add slashes to directories
if ((~isempty(data_dir)) && (data_dir(end) ~= '/'))
    data_dir = [ data_dir '/' ];
end

% define some variables which depend on the input arguments
file_ind_max = length(fnames);

% initialize variables used in the loop
filenames_all_max = 0;
filenames_all = cell(file_ind_max,1);

% loop over all matching files
for (file_ind=1:file_ind_max)
    %         % skip single frames created using the spec macro ct
    %        if (length(fnames(file_ind).name) > 7)
    %            fprintf('');
    %            if (strcmp(fnames(file_ind).name((end-6):(end-3)),'_ct.'))
    %                fprintf('skipping %s\n',fnames(file_ind).name);
    %                continue
    %            end
    %        end
    % directory: recursion
    if ((fnames(file_ind).isdir) && (recursive))
        % ignore . and .. directories
        if ((strcmp(fnames(file_ind).name,'.')) || ...
                (strcmp(fnames(file_ind).name,'..')))
            fprintf('skipping %s\n',fnames(file_ind).name);
            continue
        end
        data_dir_sub = [ data_dir fnames(file_ind).name '/' ];
        fnames_sub = dir( data_dir_sub );
        fprintf('recursion for %s\n',fnames(file_ind).name);
        [ filenames_all_rec,vararg_remain ] = ...
            collect_radial_integ_filenames(data_dir_sub, ...
            fnames_sub, ...
            integ_masks, ...
            fig_no, save_combined_I, recursive, ...
            vararg);
        % store result of this recursion
        if (~isempty(filenames_all_rec))
            filenames_all_ind = (filenames_all_max+1):(filenames_all_max+length(filenames_all_rec));
            filenames_all(filenames_all_ind) = filenames_all_rec;
            filenames_all_max = filenames_all_ind(end);
        end
        continue;
    end
    if ((length(fnames(file_ind).name) <= 4) || ...
            (strcmp(fnames(file_ind).name(end-3:end),'.tmp')) || ...
            (strcmp(fnames(file_ind).name(end-3:end),'.log')))
        fprintf('skipping %s\n',fnames(file_ind).name);
        continue
    end
    
    
    
    % store matching filenames in one array
    filenames_all_max = filenames_all_max +1;
    filenames_all{filenames_all_max} = [ data_dir fnames(file_ind).name ];
    
    
end


if  (~exist('filenames_all','var'))
    filenames_all = [];
end


if (length(filenames_all) > filenames_all_max)
    filenames_all = filenames_all{1:filenames_all_max};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [frame_I,frame_std] = ...
    perform_radial_integ(file_ind, file_ind_max, ...
    frame, ...
    integ_masks, ind_r_max, no_of_segments, use_mex, ...
    vararg)
import beamline.prep_integ_masks
import io.image_read
import plotting.plot_radial_integ
import utils.default_parameter_value
import utils.find_files
% read the raw data frame
% fprintf('%6d /%6d: ',file_ind,file_ind_max);
% [frame] = image_read(filename, vararg);
if (isempty(frame))
    error('could not load frame %u',file_ind);
end

% get the number of frames in case of multi-frame data files like HDF5
no_of_frames = size(frame,3);

% initialize result variables
frame_I = zeros(ind_r_max,no_of_segments,no_of_frames);
frame_std = zeros(ind_r_max,no_of_segments,no_of_frames);

if use_mex
    for (ind_frame = 1:no_of_frames)
        % get the current frame
        frame_data = double(frame(:,:,ind_frame));
        try
            [frame_I(:,:,ind_frame), frame_std(:,:,ind_frame)] = radial_integ_mex(int32(ind_r_max),int32(no_of_segments), integ_masks.norm_sum, integ_masks.indices, frame_data);
        catch
            tmpPath = fileparts(mfilename('fullpath'));
            fprintf('Recompiling mex function...\n');
            
            % Fall back to single thread if the OpenMP fail.
            
            eval(['mex ' fullfile(tmpPath, 'private', 'radial_integ_mex.cpp') ' -outdir ' fullfile(tmpPath, 'private')]);
            try
                [frame_I(:,:,ind_frame), frame_std(:,:,ind_frame)] = radial_integ_mex(int32(ind_r_max),int32(no_of_segments), integ_masks.norm_sum, integ_masks.indices, frame_data);
            catch ME
                fprintf('radial_integ_mex failed. If the problem persists, consider setting use_mex=false.\n');
                rethrow(ME);
            end
        end
    end
else
    for (ind_frame = 1:no_of_frames)
        
        % get the current frame
        frame_data = double(frame(:,:,ind_frame));
        %     initialize output variables for current data
        frame_I_one_frame = zeros(ind_r_max,no_of_segments);
        frame_std_one_frame = zeros(ind_r_max,no_of_segments);
        for (ind_r = 1:ind_r_max)
            for (ind_seg = 1:no_of_segments)
                if (integ_masks.norm_sum(ind_r,ind_seg) > 0)
                    frame_I_one_frame(ind_r,ind_seg) = ...
                        mean(frame_data(integ_masks.indices{ind_r,ind_seg}));
                    frame_std_one_frame(ind_r,ind_seg) = ...
                        std(frame_data(integ_masks.indices{ind_r,ind_seg}));
                else
                    % mark unknown intensities
                    frame_I_one_frame(ind_r,ind_seg) = -1;
                    frame_std_one_frame(ind_r,ind_seg) = -1;
                end
            end
        end
        frame_I(:,:,ind_frame) = frame_I_one_frame;
        frame_std(:,:,ind_frame) = frame_std_one_frame;
    end
end


%%%%%%%%%%
function [frame_I,frame_std] = ...
    perform_radial_integ_parallel(file_ind, file_ind_max, ...
    filename, ...
    integ_masks, ind_r_max, no_of_segments, ...
    vararg)
import beamline.prep_integ_masks
import io.image_read
import plotting.plot_radial_integ
import utils.default_parameter_value
import utils.find_files
% read the raw data frame
fprintf('%6d /%6d: ',file_ind,file_ind_max);
[frame] = image_read(filename, vararg);
if (isempty(frame.data))
    error('could not load %s',filename);
end

% get the number of frames in case of multi-frame data files like HDF5
no_of_frames = size(frame.data,3);

% initialize result variables
frame_I = zeros(ind_r_max,no_of_segments,no_of_frames);
frame_std = zeros(ind_r_max,no_of_segments,no_of_frames);


parfor (ind_frame = 1:no_of_frames)
    % get the current frame
    frame_data = double(frame.data(:,:,ind_frame));
    % initialize output variables for current data
    frame_I_one_frame = zeros(ind_r_max,no_of_segments);
    frame_std_one_frame = zeros(ind_r_max,no_of_segments);
    
    for (ind_r = 1:ind_r_max)
        for (ind_seg = 1:no_of_segments)
            if (integ_masks.norm_sum(ind_r,ind_seg) > 0)
                frame_I_one_frame(ind_r,ind_seg) = ...
                    mean(frame_data(integ_masks.indices{ind_r,ind_seg}));
                frame_std_one_frame(ind_r,ind_seg) = ...
                    std(frame_data(integ_masks.indices{ind_r,ind_seg}));
            else
                % mark unknown intensities
                frame_I_one_frame(ind_r,ind_seg) = -1;
                frame_std_one_frame(ind_r,ind_seg) = -1;
            end
        end
    end
    
    frame_I(:,:,ind_frame) = frame_I_one_frame;
    frame_std(:,:,ind_frame) = frame_std_one_frame;
end
