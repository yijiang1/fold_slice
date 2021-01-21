%INITIAL_CHECKS set default values for the most common variables
%
% ** p          p structure
%
% returns:
% ++ p          p structure
%
% see also: core.initialize_ptycho
%

function [p] = initial_checks(p)


import utils.*
import io.*

%%%%%%%%%%%%%
%% General %%
%%%%%%%%%%%%%

% check matlab version
check_matlab_version(9.3);

if ~isfield(p, 'use_display') || isempty(p.use_display)
    if verbose > 1
        p.use_display = true;
    else
        p.use_display = false;
    end
end

if ~usejava('desktop')
   % test if matlab was called with -nodisplay option, if yes then
   % use_display should be false
   p.use_display = false;
end

% check if fourier ptycho recon is needed
if ~isfield(p, 'fourier_ptycho')
    p.fourier_ptycho = false;
end

if ~isfield(p, 'sample_rotation_angles')
    p.sample_rotation_angles = [0,0,0];   % 3x1 vector rotation around [X,Y,beam] axes in degrees , apply a correction accounting for tilted plane oR the sample and ewald sphere curvature (high NA correction)  
end


%%% Derived quantities %%%
assert(~isempty(p.energy), 'Provide p.energy or source of metadata p.src_metadata')

% modified by YJ for electron pty
if isfield(p,'beam_source') && strcmp(p.beam_source, 'electron')
    %use relativistic corrected formula for electron pty
    p.lambda = 12.3986/sqrt((2*511.0+p.energy).*p.energy); %angstrom
else
    p.lambda = 1.23984193e-9/p.energy;                          % wavelength
end
if isscalar(p.asize); p.asize = [p.asize p.asize]; end

%%%%%%%%%%%%%%%%%%%%
%% Scan meta data %%
%%%%%%%%%%%%%%%%%%%%

% calculate fourier ptycho geometry
if p.fourier_ptycho
    if ~isfield(p, 'FP_focal_distance')
        error('For running Fourier Ptychography, please specify the focal length of your objective lens (p.FP_focal_distance)');
    end
    if ~get_option(p, 'z_lens')         
        p.z_lens = 1/(1/(p.FP_focal_distance)-1/(p.z));     
    end
end

%%%%%%%%%%%%%%%%
%% Scan queue %%
%%%%%%%%%%%%%%%%

% number of attempts to reconstruct the given dataset
if ~isfield(p.queue, 'max_attempts')
    p.queue.max_attempts = 5;
end

% lock files
if ~isfield(p.queue, 'lockfile') || isempty(p.queue.lockfile)
    if verbose > 2
        p.queue.lockfile = false;
    else
        p.queue.lockfile = true;
    end
end

if ~isfield(p.queue, 'file_queue_timeout')
    p.queue.file_queue_timeout = 10; % time to wait for a new dataset in queue
end

%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
%%%%%%%%%%%%%%%%%%%%%%

% data prefix
if isempty(p.detector.data_prefix)
    import beamline.identify_eaccount  %% not included in the ptychoshelves package 
    eaccount = identify_eaccount; 
    if ~isempty(eaccount) && eaccount(1) == 'e'
        % default setting for cSAXS beamline
        p.detector.data_prefix = [eaccount '_1_'];
    else
        verbose(0,'p.detector.data_prefix is not defined')
    end
end

% suffix for prepared data file
if ~isfield(p.prepare, 'prep_data_suffix')
    p.prepare.prep_data_suffix = '';
end

if p.asize(1) ~= p.asize(2) && (~isfield(p.prepare, 'data_preparator') || any(strcmpi(p.prepare.data_preparator, {'python', 'libDetXR','json'})))
    p.prepare.data_preparator = 'matlab_ps';
    verbose(1, 'Python preparator does not support asymmetric probe dimensions, switching to matlab_ps')
end

if p.asize(1) ~= p.asize(2) && p.prepare.force_preparation_data == false
    verbose(1, 'Loading prepared data is not supported for asymmetric p.asize, enforce load from raw data ')
    p.prepare.force_preparation_data = true; 
end


% data preparator
if ~isfield(p.prepare, 'data_preparator') || any(strcmpi(p.prepare.data_preparator, {'python', 'libDetXR','json'}))
    p.prepare.data_preparator = 'libDetXR';
    verbose(3, 'Using python data preparator.')
elseif any(strcmpi(p.prepare.data_preparator, {'matlab', 'matlab_ps','mex'}))
    p.prepare.data_preparator = 'matlab_ps';
    verbose(3, 'Using matlab data preparator.')
elseif any(strcmpi(p.prepare.data_preparator, {'matlab_aps'})) %% adde by YJ
    p.prepare.data_preparator = 'matlab_aps';
    verbose(3, 'Using matlab APS data preparator.')
else
    error('Unknown data preparator %s', p.prepare.data_preparator);
end


if ~isfield(p.prepare,'data_preparator') || isempty(p.prepare.data_preparator)
    error(' p.prepare.data_preparator is not set')
end

if strcmpi(p.prepare.data_preparator, 'matlab')
    p.prepare.data_preparator = 'matlab_ps';
elseif strcmpi(p.prepare.data_preparator, 'python')
    p.prepare.data_preparator = 'libDetXR';
end

% binning is only supported by Matlab data preparation
if isfield(p.detector,'binning')&& p.detector.binning
    p.prepare.data_preparator = 'matlab_ps';
    verbose(1, 'Using binning %ix%i, switching to matlab data loading', 2^p.detector.binning, 2^p.detector.binning)
else
    p.detector.binning = false; 
end

% binning is only supported by Matlab data preparation
if isfield(p.detector,'upsampling')&& p.detector.upsampling
    %p.prepare.data_preparator = 'matlab_ps';
    p.prepare.data_preparator = 'matlab_aps'; %modified by YJ for APS data
    verbose(1, 'Using data upsampling %ix%i, switching to matlab data loading', 2^p.detector.upsampling, 2^p.detector.upsampling)
else
    p.detector.upsampling = false; 
end

% prealignment for Fourier Ptychography
if check_option(p, 'FP_focal_distance')
    p.   fourier_ptycho = true;                                % set to true for Fourier Ptychography
else
    p.   fourier_ptycho = false;           
end

if ~isfield(p, 'prealign_FP')
    p.prealign_FP = false;
end

% Fourier Ptycho is only supported by Matlab data preparation
if p.fourier_ptycho
    p.prepare.data_preparator = 'matlab_ps';
    verbose(1, 'Switching to matlab data loading for Fourier Ptycho.')
    
    % set prealign_data to true if not distortion correction is available
    if p.prealign_FP && ~p.prealign.prealign_data && isempty(p.prealign.distortion_corr)
        p.prealign.prealign_data = true;
    end
end

% set defaults for matlab_ps
if strcmpi(p.prepare.data_preparator, 'matlab_ps')
    if ~isfield(p.io, 'data_precision')
        p.io.data_precision = 'single';
    end
    if ~isfield(p.io, 'data_nthreads')
        p.io.data_nthreads = 2;
    end
end

if isfield(p, 'prop_regime') && ~ismember(p.prop_regime, {'nearfield', 'farfield'})
    error(['Nonexistent propagation regime ', p.prop_regime ])
end    

% store prepared data
if ~isfield(p.prepare, 'store_prepared_data')
    p.prepare.store_prepared_data = true; 
end

%%%%%%%%%%%%%%%%%%%%
%% Scan positions %%
%%%%%%%%%%%%%%%%%%%%

% load positions from prepared file
if ~isfield(p.io, 'load_prep_pos')
    p.io.load_prep_pos = false;
end

%%%%%%%%%
%% I/O %%
%%%%%%%%%

% file compression
if ~isfield(p.io, 'file_compression')
    p.io.file_compression = 0;
end
if ~isfield(p.io, 'data_compression')
    p.io.data_compression = 3;
end


% run name
if ~check_option(p, 'run_name')
    % check if prefix is defined 
    if isempty(p.prefix) 
        if iscell(p.scan_str)
            p.prefix =  p.scan_str{1};
        else
            p.prefix = p.scan_str;
        end
    end
    p.run_name = sprintf('%s_%s', p.prefix, datestr(now, 'yyyy_mm_dd'));
end
verbose(3, 'run_name = %s', p.run_name);

%%%%%%%%%%%%%%%%%%%%
%% Reconstruction %%
%%%%%%%%%%%%%%%%%%%%

% backward compatibilty for initial_iterate
if isfield(p, 'initial_iterate') && ~isfield(p, 'initial_iterate_object')
    p.initial_iterate_object = p.initial_iterate;
    p = rmfield(p, 'initial_iterate');
end
if isfield(p, 'initial_iterate_file') && ~isfield(p, 'initial_iterate_object_file')
    p.initial_iterate_object_file = p.initial_iterate_file;
    p = rmfield(p, 'initial_iterate_file');
end

% model probe
if ~isfield(p.model, 'probe_central_stop')
    p.model.probe_central_stop = false;
end

if ~isfield(p.model, 'probe_central_stop_diameter') && p.model.probe_central_stop
    p.model.probe_central_stop_diameter = 50e-6;
end

%%%%%%%%%%%%%%%%%%%
%% Plot and save %%
%%%%%%%%%%%%%%%%%%%

% plot prepared data
if ~isfield(p.plot, 'prepared_data') || (isfield(p.plot, 'prepared_data')&& isempty(p.plot.prepared_data))
    if p.verbose_level > 2
        p.plot.prepared_data = true;
    else
        p.plot.prepared_data = false;
    end
end

% plotting
if ~isfield(p.plot, 'interval') || isempty(p.plot.interval)
    if verbose > 2
        p.plot.interval = 10;
    else
        p.plot.interval = 200;
    end
end

% external call to save figures
if ~isfield(p.save, 'external')
    p.save.external = false;
end

% propagation and apodization
if ~isfield(p.plot, 'obj_apod')
    p.plot.obj_apod = false;
end
if ~isfield(p.plot, 'prop_obj')
    p.plot.prop_obj = 0;
end

% calculate FSC
if ~isfield(p.plot, 'calc_FSC')
    p.plot.calc_FSC = false;
end
if ~isfield(p.plot, 'show_FSC')
    p.plot.show_FSC = utils.verbose>2;
end
if ~isfield(p.plot, 'probe_spectrum')|| isempty(p.plot.probe_spectrum)
    p.plot.probe_spectrum = utils.verbose>2;
end
if ~isfield(p.plot, 'object_spectrum')|| isempty(p.plot.object_spectrum)
    p.plot.object_spectrum = utils.verbose>2;
end
if ~isfield(p.save, 'store_images_ids' )|| isempty(p.save.store_images_ids)
    p.save.store_images_ids = 1:4; 
end


%%%%%%%%%%%%%
%% Engines %%
%%%%%%%%%%%%%

% at least one engine has to be specified
if ~isfield(p, 'engines')
    error('At least one reconstruction engine has to be selected. Please check your template.')
end



% first engine is external
p.external_engine0 = strcmpi(p.engines{1}.name, 'c_solver') || ...
    (isfield(p.engines{1}, 'external') && p.engines{1}.external);

if ~isfield(p,'remove_scaling_ambiguity')
    % if true. try to keep norm(probe) constant during the reconstruction
    p.remove_scaling_ambiguity = true; 
end


end

