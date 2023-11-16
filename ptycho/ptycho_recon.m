% Wrapper function for reconstruction using the GPU engines
% Can be called in python and used for automatic parameter tuning
% Written by Yi Jiang
%
% Input:
% ** param    - structure containting common parameters for reconstruction

% Outputs:
% ++ out           - reconstruction
% ++ eng           - engine structure
% ++ recon_error   - averaged data error of last iteration

function [out, eng, data_error] = ptycho_recon(param)

    addpath(strcat(pwd,'/utils/'))
    addpath(core.find_base_package)

    % parse inputs
    parse_param = inputParser;
    parse_param.KeepUnmatched = true;
    
    % basic recon parameters
    parse_param.addParameter('eng_name', 'GPU', @ischar)
    parse_param.addParameter('Niter', 1000, @isnumeric)
    parse_param.addParameter('Niter_save_results_every', 1000, @isnumeric)
    parse_param.addParameter('Nprobe', 1, @isnumeric)
    parse_param.addParameter('grouping', 100, @isnumeric)
    parse_param.addParameter('time_limit', inf, @isnumeric)
    parse_param.addParameter('gpu_id', 1, @isnumeric)
    parse_param.addParameter('background', 0, @isnumeric)
    parse_param.addParameter('apply_multimodal_update', 0, @islogical)
    parse_param.addParameter('probe_position_search', inf, @isnumeric)
    parse_param.addParameter('apply_geometry_model_constraint', 0, @islogical)
    parse_param.addParameter('geometry_model_scale', 0, @islogical)
    parse_param.addParameter('geometry_model_shear', 0, @islogical)
    parse_param.addParameter('geometry_model_rotation', 0, @islogical)
    parse_param.addParameter('geometry_model_asymmetry', 0, @islogical)
    parse_param.addParameter('update_pos_weight_every', inf, @isnumeric)

    parse_param.addParameter('variable_probe_modes', 0, @isnumeric)
    parse_param.addParameter('variable_probe_smooth', 0, @isnumeric)
    parse_param.addParameter('variable_intensity', 0, @islogical)
    parse_param.addParameter('method',  'MLs', @ischar)
    parse_param.addParameter('momentum',  0, @isnumeric)
    parse_param.addParameter('delta_p',  0.1, @isnumeric)
    parse_param.addParameter('beta_probe',  1, @isnumeric)
    parse_param.addParameter('beta_object',  1, @isnumeric)
    parse_param.addParameter('opt_errmetric',  'L1', @ischar)
    parse_param.addParameter('reg_mu',  0, @isnumeric)
    parse_param.addParameter('positivity_constraint_object', 0, @isnumeric)
    parse_param.addParameter('Ndp_presolve', inf, @isnumeric)
    parse_param.addParameter('diff_pattern_blur', 0, @isnumeric)
    
    parse_param.addParameter('dp_custom_fliplr', 0, @isnumeric)
    parse_param.addParameter('dp_custom_flipud', 0, @isnumeric)
    parse_param.addParameter('dp_custom_transpose', 0, @isnumeric)
    
    %multislice ptycho parameters
    parse_param.addParameter('Nlayers', 2, @isnumeric)
    parse_param.addParameter('delta_z', 10, @isnumeric)
    parse_param.addParameter('regularize_layers', 0.1, @isnumeric)
    parse_param.addParameter('init_layer_preprocess', '', @ischar)
    parse_param.addParameter('init_layer_append_mode', '', @ischar)
    parse_param.addParameter('init_layer_scaling_factor', 1, @isnumeric)
    parse_param.addParameter('layer4pos', 0, @isnumeric)

    % data I/O parameters
    parse_param.addParameter('base_path', '', @ischar)
    parse_param.addParameter('scan_number', 1 , @isnumeric)
    parse_param.addParameter('beam_source', '', @ischar) %for e-ptycho.

    parse_param.addParameter('Ndp', 256, @isnumeric )
    parse_param.addParameter('cen_dp_y', 0, @isnumeric)
    parse_param.addParameter('cen_dp_x', 0, @isnumeric)

    parse_param.addParameter('avg_photon_threshold', 0, @isnumeric)

    parse_param.addParameter('detector_dist', 1.0, @isnumeric)
    parse_param.addParameter('energy', 8.8, @isnumeric)
    parse_param.addParameter('d_alpha', 1.0, @isnumeric ) %for e-ptycho.
    parse_param.addParameter('dk', 0, @isnumeric ) %for e-ptycho.. Note: dk is checked before d_alpha

    parse_param.addParameter('affine_matrix_11', 1, @isnumeric)
    parse_param.addParameter('affine_matrix_12', 0, @isnumeric)
    parse_param.addParameter('affine_matrix_21', 0, @isnumeric)
    parse_param.addParameter('affine_matrix_22', 1, @isnumeric)

    parse_param.addParameter('angular_correction_setup', 'none', @ischar)

    parse_param.addParameter('use_model_object', 1, @islogical)
    parse_param.addParameter('initial_object_file', '', @ischar)
    parse_param.addParameter('multiple_layers_obj', 0, @islogical)
    parse_param.addParameter('sum_obj_layers', 0, @islogical)

    parse_param.addParameter('use_model_probe', 0, @islogical)
    parse_param.addParameter('initial_probe_file', '', @ischar)
    parse_param.addParameter('normalize_init_probe', 0, @islogical)
    parse_param.addParameter('crop_pad_init_probe', 0, @islogical)
    parse_param.addParameter('probe_file_propagation', 0, @isnumeric)
    parse_param.addParameter('manual_center_probe_x', 0, @isnumeric)
    parse_param.addParameter('manual_center_probe_y', 0, @isnumeric)

    parse_param.addParameter('model_probe_prop_dist', 0, @isnumeric)
    parse_param.addParameter('model_probe_outer_zone_width', 60e-9, @isnumeric)
    parse_param.addParameter('model_probe_zone_plate_diameter', 114.8e-6, @isnumeric)
    
    % for model STEM probe
    parse_param.addParameter('probe_alpha_max', 20, @isnumeric)
    parse_param.addParameter('probe_df', 0, @isnumeric)
    parse_param.addParameter('probe_c3', 0, @isnumeric)
    parse_param.addParameter('probe_c5', 0, @isnumeric)
    parse_param.addParameter('probe_c7', 0, @isnumeric)
    parse_param.addParameter('probe_f_a2', 0, @isnumeric)
    parse_param.addParameter('probe_theta_a2', 0, @isnumeric)
    parse_param.addParameter('probe_f_a3', 0, @isnumeric)
    parse_param.addParameter('probe_theta_a3', 0, @isnumeric)
    parse_param.addParameter('probe_f_c3', 0, @isnumeric)
    parse_param.addParameter('probe_theta_c3', 0, @isnumeric)
    parse_param.addParameter('scan_format', '', @ischar)
    parse_param.addParameter('scan_string_format', '', @ischar)
    parse_param.addParameter('roi_label', '', @ischar)
    
    parse_param.addParameter('src_metadata', 'none', @ischar)
    parse_param.addParameter('detector_name', 'eiger_APS', @ischar)
    parse_param.addParameter('detector_upsampling', 0, @isnumeric)
    
    parse_param.addParameter('data_preparator', 'matlab_aps', @ischar)
    parse_param.addParameter('src_positions', 'hdf5_pos', @ischar)
    parse_param.addParameter('positions_file', '', @ischar)
    parse_param.addParameter('scan_type', 'default', @ischar)
    parse_param.addParameter('custom_positions_source', '', @ischar)
    
    parse_param.addParameter('scan_nx', 1, @isnumeric)
    parse_param.addParameter('scan_ny', 1, @isnumeric)
    parse_param.addParameter('scan_step_size_x', 1, @isnumeric)
    parse_param.addParameter('scan_step_size_y', 1, @isnumeric)
    parse_param.addParameter('scan_custom_fliplr', 0, @isnumeric)
    parse_param.addParameter('scan_custom_flipud', 0, @isnumeric)
    parse_param.addParameter('scan_custom_transpose', 0, @isnumeric)
    
    parse_param.addParameter('output_dir_base','', @ischar)
    parse_param.addParameter('output_dir_prefix','', @ischar)
    parse_param.addParameter('output_dir_suffix','', @ischar)
    
    parse_param.addParameter('save_init_probe', 0, @islogical)

    parse_param.addParameter('rng_seed', -1, @isnumeric)
    parse_param.addParameter('verbose_level',  2 , @isnumeric)
    
    parse_param.parse(param)
    param_input = parse_param.Results;
    
    %% check inputs
    assert(~isempty(param_input.base_path), 'Base path cannot be empty!')
    assert(~isempty(param_input.output_dir_base), 'Output path cannot be empty!')

    assert(param_input.Niter>=1, 'Invalid number of iterations!')
    assert(param_input.Nprobe>=1, 'Invalid number of probes!')
    assert(param_input.grouping>1, 'Invalid group size!')

    assert(param_input.probe_alpha_max>0, 'Invalid aperature size!')
    
    %% set rng seed
    if param_input.rng_seed>=0
        rng(param_input.rng_seed);
    end
    %%
    Nfly = 0;
    Ndp = double(param_input.Ndp);  % size of cbed
    if param_input.cen_dp_y>0
        cen_dp_y = param_input.cen_dp_y;
    else
        cen_dp_y = floor(Ndp/2) + 1;
    end
    if param_input.cen_dp_x>0
        cen_dp_x = param_input.cen_dp_x;
    else
        cen_dp_x = floor(Ndp/2) + 1;
    end
    
    %% General
    p = struct();
    p.   verbose_level = param_input.verbose_level;                            % verbosity for standard output (0-1 for loops, 2-3 for testing and adjustments, >= 4 for debugging)
    p.   use_display = false;                                      % global switch for display, if [] then true for verbose > 1
    p.   scan_number = param_input.scan_number;                                    % Multiple scan numbers for shared scans

    p.   avg_photon_threshold = param_input.avg_photon_threshold;

    % Geometry
    p.   z = param_input.detector_dist;                         % Distance from object to detector 
    p.   asize = [Ndp, Ndp];                                    % Diffr. patt. array size   
    p.   ctr = [cen_dp_y, cen_dp_x];                                % Diffr. patt. center coordinates (y,x) (empty means middle of the array); e.g. [100 207;100+20 207+10];
    p.   beam_source = param_input.beam_source;                 % Added by YJ for electron pty. Use relativistic corrected formula for wavelength. Also change the units on figures
    p.   d_alpha = param_input.d_alpha;                        % Added by YJ. d_alpha is the pixel size in cbed (mrad). This is used to determine pixel size in electron ptycho
    p.   dk = param_input.dk;                                   % Added by YJ. dk is the pixel size in cbed (1/ansgtrom). This is used to determine pixel size in electron ptycho
    
    p.   prop_regime = 'farfield';                              % propagation regime: nearfield, farfield (default), !! nearfield is supported only by GPU engines 
    p.   focus_to_sample_distance = [];                         % sample to focus distance, parameter to be set for nearfield ptychography, otherwise it is ignored 
    p.   FP_focal_distance = [];                                %  if nonempty -> assume Fourier ptychography configuration, FP_focal_distance = focal length of objective lens for Fourier Ptychography only,
    p.   angular_correction_setup = param_input.angular_correction_setup;                         % if src_positions=='orchestra', choose angular correction for specific cSAXS experiment: 'flomni', 'omny', 'lamni', 'none', 
    if param_input.energy < 0                                            % Energy (in keV), leave empty to use spec entry mokev
        p.   energy = [];                                           % Energy (in keV), leave empty to use spec entry mokev
    else
        p.   energy = param_input.energy;                                           % Energy (in keV), leave empty to use spec entry mokev
    end
    p.   sample_rotation_angles = [0,0,0];                      % Offaxis ptychography correction , 3x1 vector rotation around [X,Y,beam] axes in degrees , apply a correction accounting for tilted plane oR the sample and ewald sphere curvature (high NA correction)

    %p.   affine_angle = 0;                                     % Not used by ptycho_recons at all. This allows you to define a variable for the affine matrix below and keep it in p for future record. This is used later by the affine_matrix_search.m script
	p.   affine_matrix = [param_input.affine_matrix_11, param_input.affine_matrix_12;...
                          param_input.affine_matrix_21, param_input.affine_matrix_22];
    % Scan meta data
    p.   src_metadata = param_input.src_metadata;                                 % source of the meta data, following options are supported: 'spec', 'none' , 'artificial' - or add new to +scan/+meta/

    % Scan queue
    p.   queue.name = '';                                       % specify file queue; currently only 'filelist' is supported
    p.   queue.path=['//'];      % Folder where the queue of files is defined, note the content of files can overwrite some parameters in p-structure
    p.   queue.max_attempts = 5;                                % Max number of attempts to reconstruct a scan.
    p.   queue.file_queue_timeout = 10;                         % Time to wait when queue is empty before checking it again 
    p.   queue.remote_recons = false;                           % divide the reconstruction into primary/replica processes to reconstruction on a remote server
    p.   queue.recon_latest_first = 1;                          % When using 'p.queue_path', (1) reconstruct the latest measurement first or (0) reconstruct in lexicographical order
    p.   queue.remote_path = '';                                % Queue list for remote reconstructions. Needs to be accessible for primary and replica processes
    p.   queue.tmp_dir_remote = '';                             % shared directory for storing the remote reconstruction

    p.   queue.lockfile = false;                                % If true writes a lock file, if lock file exists skips recontruction
    p.   spec.waitforscanfinish = true;                         % Checks spec file for the scan end flag 'X#'
    p.   spec.check_nextscan_started = true;                    % Waits until the next scan starts to begin reconstructing this one. It is important for OMNY scans with orchestra
    p.   spec.isptycho = {};                                    % Use only when SPEC is used: = {'round_roi','cont_line','ura_mesh'}  ( = {} to skip)  List of ptycho spec commands for valid ptycho scans

    % Data preparation
    p.   detector.name = param_input.detector_name;                           % see +detectors/ folder 
    p.   detector.check_2_detpos = [];                          % = []; (ignores)   = 270; compares to dettrx to see if p.ctr should be reversed (for OMNY shared scans 1221122), make equal to the middle point of dettrx between the 2 detector positions
    p.   detector.data_prefix = '';                             % Default using current eaccount e.g. e14169_1_
    p.   detector.binning = false;                              % = true to perform 2x2 binning of detector pixels, for binning = N do 2^Nx2^N binning
    p.   detector.upsampling = param_input.detector_upsampling;  % upsample the measured data by 2^data_upsampling, (transposed operator to the binning), it can be used for superresolution in nearfield ptychography or to account for undersampling in a far-field dataset
    p.   detector.burst_frames = 1;                             % number of frames collected per scan position

    p.   prepare.data_preparator = param_input.data_preparator;                    % data preparator; 'python' or 'matlab' 
    p.   prepare.auto_prepare_data = true;                      % if true: prepare dataset from raw measurements if the prepared data does not exist
    p.   prepare.force_preparation_data = true;                 % Prepare dataset even if it exists, it will overwrite the file % Default: @prepare_data_2d
    p.   prepare.store_prepared_data = false;                    % store the loaded data to h5 even for non-external engines (i.e. other than c_solver)
    p.   prepare.prepare_data_function = '';                    % (used only if data should be prepared) custom data preparation function handle;
    p.   prepare.auto_center_data = false;                      % if matlab data preparator is used, try to automatically center the diffraction pattern to keep center of mass in center of diffraction

    p.   prealign_FP = false;                                   % use prealignment routines for Fourier Ptychography
    p.   prealign.asize = [1000 1000];                          % array size for the alignment procedure
    p.   prealign.crop_dft = 100;                               % crop the dftregistration input
    p.   prealign.prealign_data = true;                         % recalculate the alignment
    p.   prealign.axis = 1;                                     % alignment axis
    p.   prealign.type = {'round'};                             % alignment routine
    p.   prealign.numiter = 5;                                  % number of alignment iterations
    p.   prealign.rad_filt_min = 25e-6;                         % discard positions < rad_filt_min radius
    p.   prealign.rad_filt_max = 80e-6;                         % discard positions > rad_filt_max radius
    p.   prealign.load_alignment = true;                        % load alignment from an alignment_file
    p.   prealign.alignment_file = '';      % alignment file
    p.   prealign.mag_est = 160;                                % estimated magnification; used as an initial guess for the distortion correction matrix
    p.   prealign.use_distortion_corr = true;                   % use distortion correction; if distortion_corr is empty, it will calculate a new correction based on the shifts retrieved from the alignment
    p.   prealign.distortion_corr = [];                         % distortion correction matrix; [161.3003, 3.4321, -6.7294, 0.0000, 0.9675, 2.0220, 0.0540];

    % Scan positions
    p.   src_positions = param_input.src_positions;                           % 'spec', 'orchestra', 'load_from_file', 'matlab_pos' (scan params are defined below) or add new position loaders to +scan/+positions/
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p.   positions_file = [param_input.positions_file];    % Filename pattern for position files, Example: ['../../specES1/scan_positions/scan_%05d.dat']; (the scan number will be automatically filled in)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p.   spec.motor.fine_motors = {};                           % Y and X motor name for positions, leave empty for defaults
    p.   spec.motor.fine_motors_scale = [];                     % ptycho expects real positions in m; 
    p.   spec.motor.coarse_motors = {};                         % Coarse sample position for shared object, use {X-motor, Y-motor} 
    p.   spec.motor.coarse_motors_scale = [];                   % Scale of the coarse motors (to scale the provided values to meters)

    % scan parameters for option src_positions = 'matlab_pos';
    p.   scan.type = param_input.scan_type;                               % {'round', 'raster', 'round_roi', 'custom'}
    p.   scan.roi_label = param_input.roi_label;                            % For APS data
    p.   scan.format = param_input.scan_format;                      % For APS data format for scan directory generation

    p.   scan.radius_in = 0;                                    % round scan: interior radius of the round scan
    p.   scan.radius_out = 5e-6;                                % round scan: exterior radius of the round scan
    p.   scan.nr = 10;                                          % round scan: number of intervals (# of shells - 1)
    p.   scan.nth = 3;                                          % round scan: number of points in the first shell
    p.   scan.lx = 20e-6;                                       % round_roi scan: width of the roi
    p.   scan.ly = 20e-6;                                       % round_roi scan: height of the roi
    p.   scan.dr = 1.5e-6;                                      % round_roi scan: shell step size
    p.   scan.nx = param_input.scan_nx;                         % raster scan: number of steps in x
    p.   scan.ny = param_input.scan_ny;                         % raster scan: number of steps in y
    p.   scan.step_size_x = param_input.scan_step_size_x;       % raster scan: step size (grid spacing)
    p.   scan.step_size_y = param_input.scan_step_size_y;       % raster scan: step size (grid spacing)
    p.   scan.custom_flip = [param_input.scan_custom_fliplr,...
                             param_input.scan_custom_flipud,...
                             param_input.scan_custom_transpose]; % raster scan: apply custom flip [fliplr, flipud, transpose] to positions- similar to eng.custom_data_flip in GPU engines. Added by ZC.

    p.   scan.step_randn_offset = 0;                            % raster scan: relative random offset from the ideal periodic grid to avoid the raster grid pathology 
    p.   scan.b = 0;                                            % fermat: angular offset
    p.   scan.n_max = 1e4;                                      % fermat: maximal number of points generated 
    p.   scan.step = 0.5e-6;                                      % fermat: step size 
    p.   scan.cenxy = [0,0];                                    % fermat: position of center offset 
    p.   scan.roi = [];                                         % Region of interest in the object [xmin xmax ymin ymax] in meters. Points outside this region are not used for reconstruction.
                                                                %  (relative to upper corner for raster scans and to center for round scans)   
    p.   scan.custom_positions_source = param_input.custom_positions_source;                     % custom: a string name of a function that defines the positions; also accepts mat file with entry 'pos', see +scans/+positions/+mat_pos.m
    p.   scan.custom_params = [];                               % custom: the parameters to feed to the custom position function.

    % I/O
    p.   prefix = '';                                           % For automatic output filenames. If empty: scan number
    p.   suffix = '';                                      % Optional suffix for reconstruction 
    p.   scan_string_format = param_input.scan_string_format;                          % format for scan string generation, it is used e.g for plotting and data saving 
    p.   base_path = strcat(param_input.base_path);                                  % base path : used for automatic generation of other paths 

    p.   specfile = '';                                         % Name of spec file to get motor positions and check end of scan, defaut is p.spec_file == p.base_path;
    p.   ptycho_matlab_path = pwd;                               % cSAXS ptycho package path
    p.   cSAXS_matlab_path = '';                                % cSAXS base package path
    p.   ptycho_package_path = pwd;                              % ptycho package path. If empty it uses fullfile(p.base_path,

    p.   raw_data_path{1} = '';                                 % Default using compile_x12sa_filename, used only if data should be prepared automatically
    p.   prepare_data_path = '';                                % Default: base_path + 'analysis'. Other example: '/afs/psi.ch/project/CDI/cSAXS_project/analysis2/'; also supports %u to insert the scan number at a later point (e.g. '/afs/psi.ch/project/CDI/cSAXS_project/analysis2/S%.5u')
    p.   prepare_data_filename = [];                            % Leave empty for default file name generation, otherwise use [sprintf('S%05d_data_%03dx%03d',p.scan_number(1), p.asize(1), p.asize(2)) p.prep_data_suffix '.h5'] as default 
    p.   save_path{1} = '';                                     % Default: base_path + 'analysis'. Other example: '/afs/psi.ch/project/CDI/cSAXS_project/analysis2/'; also supports %u to insert the scan number at a later point (e.g. '/afs/psi.ch/project/CDI/cSAXS_project/analysis2/S%.5u')
    p.   io.default_mask_file = '';                             % load detector mask defined in this file instead of the mask in the detector packages, (used only if data should be prepared) 
    p.   io.default_mask_type = 'binary';                       % (used only if data should be prepared) ['binary', 'indices']. Default: 'binary' 
    p.   io.file_compression = 0;                               % reconstruction file compression for HDF5 files; 0 for no compression
    p.   io.data_compression = 3;                               % prepared data file compression for HDF5 files; 0 for no compression
    p.   io.load_prep_pos = false;                              % load positions from prepared data file and ignore positions provided by metadata

    p.   io.data_descriptor = '';                    % added by YJ. A short string that describe data when sending notifications 
    p.   io.phone_number = [];                      % phone number for sending messages
    p.   io.send_failed_scans_SMS = false;                      % send message if p.queue_max_attempts is exceeded
    p.   io.send_finished_recon_SMS = false;                     % send message after the reconstruction is completed
    p.   io.send_crashed_recon_SMS = false;                     % send message if the reconstruction crashes
    p.   io.SMS_sleep = 1800;                                   % max 1 message per SMS_sleep seconds
    p.   io.script_name = mfilename;                            % added by YJ. store matlab script name

    p.   artificial_data_file = 'template_artificial_data';     % artificial data parameters, set p.src_metadata = 'artificial' to use this template

    %% Reconstruction
    % Initial iterate object
    p.   model_object = param_input.use_model_object;           % Use model object, if false load it from file 
    p.   model.object_type = 'rand';                            % specify how the object shall be created; use 'rand' for a random initial guess; use 'amplitude' for an initial guess based on the prepared data
    p.   multiple_layers_obj = param_input.multiple_layers_obj;   % Keep all object layers from a multislice reconstruction
    p.   sum_obj_layers = param_input.sum_obj_layers;             % Sum all object layers from a multislice reconstruction
    p.   initial_iterate_object_file{1} = param_input.initial_object_file;  %  use this mat-file as initial guess of object, it is possible to use wild characters and pattern filling, example: '../analysis/S%05i/wrap_*_1024x1024_1_recons*'

    % Initial iterate probe
    p.   model_probe = param_input.use_model_probe;                                    % Use model probe, if false load it from file 
    p.   model.probe_is_focused = true;                         % Model probe is focused (false: just a pinhole)
    p.   model.probe_central_stop = true;                       % Model central stop
    p.   model.probe_diameter = 170e-6;                         % Model probe pupil diameter
    p.   model.probe_central_stop_diameter = 50e-6;             % Model central stop diameter
    p.   model.probe_zone_plate_diameter = param_input.model_probe_zone_plate_diameter;              % Model probe zone plate diameter
    p.   model.probe_outer_zone_width = param_input.model_probe_outer_zone_width;                     % Model probe zone plate outermost zone width (not used if not a focused probe) 
    p.   model.probe_propagation_dist = param_input.model_probe_prop_dist;                 % Model probe propagation distance (pinhole <-> sample for unfocused, focal-plane <-> sample for focused)
    p.   model.probe_focal_length = 51e-3;                      % Model probe focal length (used only if model_is_focused is true
                                                                %   AND model_outer_zone_width is empty)
    p.   model.probe_upsample = 10;                             % Model probe upsample factor (for focused probes)

    p.   model.probe_alpha_max = param_input.probe_alpha_max;   % Model STEM probe's aperture size
    p.   model.probe_df = param_input.probe_df;                	% Model STEM probe's defocus
    p.   model.probe_c3 = param_input.probe_c3;                 % Model STEM probe's third-order spherical aberration in angstrom (optional)
    p.   model.probe_c5 = param_input.probe_c5;                 % Model STEM probe's fifth-order spherical aberration in angstrom (optional)
    p.   model.probe_c7 = param_input.probe_c7;                 % Model STEM probe's seventh-order spherical aberration in angstrom (optional)
    p.   model.probe_f_a2 = param_input.probe_f_a2;             % Model STEM probe's twofold astigmatism in angstrom
    p.   model.probe_theta_a2 = param_input.probe_theta_a2;     % Model STEM probe's twofold azimuthal orientation in radian
    p.   model.probe_f_a3 = param_input.probe_f_a3;             % Model STEM probe's threefold astigmatism in angstrom
    p.   model.probe_theta_a3 = param_input.probe_theta_a3;     % Model STEM probe's threefold azimuthal orientation in radian
    p.   model.probe_f_c3 = param_input.probe_f_c3;             % Model STEM probe's coma in angstrom
    p.   model.probe_theta_c3 = param_input.probe_theta_c3;     % Model STEM probe's coma azimuthal orientation in radian

    p.   initial_probe_file = param_input.initial_probe_file;        % Use probe from this mat-file (not used if model_probe is true)
    p.   probe_file_propagation = param_input.probe_file_propagation;  % Distance for propagating the probe from file in meters (xray) or in angstroms (electron), = 0 to ignore
    p.   normalize_init_probe = param_input.normalize_init_probe;  % Added by YJ. Can be used to disable normalization of initial probes
    p.   crop_pad_init_probe = param_input.crop_pad_init_probe;    % added by YJ. Crop/pad initial probe if size is different. Default is false and interpolation is used. Useful when diffraction pattern is upsampled.

    % Shared scans - Currently working only for sharing probe and object
    p.   share_probe  = 0;                                       % Share probe between scans. Can be either a number/boolean or a list of numbers, specifying the probe index; e.g. [1 2 2] to share the probes between the second and third scan. 
    p.   share_object = 0;                                      % Share object between scans. Can be either a number/boolean or a list of numbers, specifying the object index; e.g. [1 2 2] to share the objects between the second and third scan. 

    % Modes
    p.   probe_modes  = param_input.Nprobe;                                      % Number of coherent modes for probe
    p.   object_modes = 1;                                      % Number of coherent modes for object
    % Mode starting guess
    p.   mode_start_pow = [0.02];                               % Normalized intensity on probe modes > 1. Can be a number (all higher modes equal) or a vector
    p.   mode_start = 'herm';                                   % (for probe) = 'rand', = 'herm' (Hermitian-like base), = 'hermver' (vertical modes only), = 'hermhor' (horizontal modes only)
    p.   ortho_probes = true;                                   % orthogonalize probes after each engine

    %% Plot, save and analyze
    p.   plot.prepared_data = false;                         % plot prepared data
    p.   plot.interval = inf;                   % plot each interval-th iteration, does not work for c_solver code
    p.   plot.log_scale = [0 0];                                % Plot on log scale for x and y
    p.   plot.realaxes = true;                                  % Plots show scale in microns
    p.   plot.remove_phase_ramp = false;                        % Remove phase ramp from the plotted / saved phase figures 
    p.   plot.fov_box = false;                                   % Plot the scanning FOV box on the object (both phase and amplitude)
    p.   plot.fov_box_color = 'r';                              % Color of the scanning FOV box
    p.   plot.positions = true;                                 % Plot the scanning positions
    p.   plot.mask_bool = true;                                 % Mask the noisy contour of the reconstructed object in plots
    p.   plot.windowautopos = true;                             % First plotting will auto position windows
    p.   plot.obj_apod = false;                                 % Apply apodization to the reconstructed object;
    p.   plot.prop_obj = 0;                                     % Distance to propagate reconstructed object before plotting [m]
    p.   plot.show_layers = true;                               % show each layer in multilayer reconstruction 
    p.   plot.show_layers_stack = false;                        % show each layer in multilayer reconstruction by imagesc3D
    p.   plot.object_spectrum = [];                             % Plot propagated object (FFT for conventional ptycho); if empty then default is false if verbose_level < 3 and true otherwise
    p.   plot.probe_spectrum = [];                              % Plot propagated probe (FFT for conventional ptycho); if empty then default is false if verbose_level < 3 and true otherwise
    p.   plot.conjugate = false;                                % plot complex conjugate of the reconstruction 
    p.   plot.horz_fact = 2.5;                                  % Scales the space that the ptycho figures take horizontally
    p.   plot.FP_maskdim = 180e-6;                              % Filter the backpropagation (Fourier Ptychography)
    p.   plot.calc_FSC = false;                                 % Calculate the Fourier Shell correlation for 2 scans or compare with model in case of artificial data tests 
    p.   plot.show_FSC = false;                                 % Show the FSC plots, including the cropped FOV
    p.   plot.residua = true;                                  % highlight phase-residua in the image of the reconstructed phase

    p.   save.external = false;                             % Use a new Matlab session to run save final figures (saves ~6s per reconstruction). Please be aware that this might lead to an accumulation of Matlab sessions if your single reconstruction is very fast.
    p.   save.store_images = false;                              % Write preview images containing the final reconstructions in [p.base_path,'analysis/online/ptycho/'] if p.use_display = 0 then the figures are opened invisible in order to create the nice layout. It writes images in analysis/online/ptycho
    p.   save.store_images_intermediate = false;                % save images to disk after each engine
    p.   save.store_images_ids = 1:4;                           % identifiers  of the figure to be stored, 1=obj. amplitude, 2=obj. phase, 3=probes, 4=errors, 5=probes spectrum, 6=object spectrum
    p.   save.store_images_format = 'png';                      % data type of the stored images jpg or png 
    p.   save.store_images_dpi = 150;                           % DPI of the stored bitmap images 
    p.   save.exclude = {'fmag', 'fmask', 'illum_sum'};         % exclude variables to reduce the file size on disk
    p.   save.save_reconstructions_intermediate = false;        % save final object and probes after each engine
    p.   save.save_reconstructions = false;                      % save reconstructed object and probe when full reconstruction is finished 
    p.   save.output_file = 'h5';                               % data type of reconstruction file; 'h5' or 'mat'

    %% ENGINES
    % --------- GPU engines  -------------   See for more details: Odstrčil M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
    eng = struct();                        % reset settings for this engine 
    eng. name = param_input.eng_name;    
    eng. use_gpu = true;                   % if false, run CPU code, but it will get very slow 
    eng. keep_on_gpu = true;               % keep data + projections on GPU, false is useful for large data if DM is used
    eng. compress_data = false;             % use automatic online memory compression to limit need of GPU memory
    eng. gpu_id = param_input.gpu_id;                      % default GPU id, [] means choosen by matlab
    eng. check_gpu_load = true;            % check available GPU memory before starting GPU engines 

    % general 
    eng. time_limit = param_input.time_limit;
    eng. number_iterations = param_input.Niter;          % number of iterations for selected method 
    if Ndp > double(param_input.Ndp_presolve)
        eng. asize_presolve = [param_input.Ndp_presolve, param_input.Ndp_presolve];      % crop data to "asize_presolve" size to get low resolution estimate that can be used in the next engine as a good initial guess 
    else
        eng. asize_presolve = [];      % crop data to "asize_presolve" size to get low resolution estimate that can be used in the next engine as a good initial guess 
    end
    eng. share_probe = p.share_probe;                 % Share probe between scans. Can be either a number/boolean or a list of numbers, specifying the probe index; e.g. [1 2 2] to share the probes between the second and third scan.
    eng. share_object = p.share_object;                % Share object between scans. Can be either a number/boolean or a list of numbers, specifying the object index; e.g. [1 2 2] to share the objects between the second and third scan. 
    eng. align_shared_objects = false;     % before merging multiple unshared objects into one shared, the object will be aligned and the probes shifted by the same distance -> use for alignement and shared reconstruction of drifting scans  

    eng. method = param_input.method;                   % choose GPU solver: DM, ePIE, hPIE, MLc, Mls, -- recommended are MLc and MLs
    eng. opt_errmetric = param_input.opt_errmetric;            % optimization likelihood - poisson, L1
    eng. grouping = param_input.grouping;                   % size of processed blocks, larger blocks need more memory but they use GPU more effeciently, !!! grouping == inf means use as large as possible to fit into memory 
                                           % * for hPIE, ePIE, MLs methods smaller blocks lead to faster convergence, 
                                           % * for MLc the convergence is similar 
                                           % * for DM is has no effect on convergence
    %eng. probe_modes  = 1;                % Number of coherent modes for probe
    eng. object_change_start = 1;          % Start updating object at this iteration number
    eng. probe_change_start = 1;           % Start updating probe at this iteration number
    eng. probe_modes  = p.probe_modes;                % Number of coherent modes for probe

    % regularizations
    eng. reg_mu = param_input.reg_mu;                       % Regularization (smooting) constant ( reg_mu = 0 for no regularization)
    eng. delta = 0;                        % press values to zero out of the illumination area in th object, usually 1e-2 is enough 
    eng. positivity_constraint_object = param_input.positivity_constraint_object; % enforce weak (relaxed) positivity in object, ie O = O*(1-a)+a*|O|, usually a=1e-2 is already enough. Useful in conbination with OPRP or probe_fourier_shift_search  

    eng. apply_multimodal_update = param_input.apply_multimodal_update; % apply all incoherent modes to object, it can cause isses if the modes collect some crap 
    eng. probe_backpropagate = 0;         % backpropagation distance the probe mask, 0 == apply in the object plane. Useful for pinhole imaging where the support can be applied  at the pinhole plane
    eng. probe_support_radius = [];       % Normalized radius of circular support, = 1 for radius touching the window    
    eng. probe_support_fft = false;       % assume that there is not illumination intensity out of the central FZP cone and enforce this contraint. Useful for imaging with focusing optics. Helps to remove issues from the gaps between detector modules.

    % basic recontruction parameters 
    % PIE / ML methods                    % See for more details: Odstrčil M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
    eng. beta_object = param_input.beta_object;	% object step size, larger == faster convergence, smaller == more robust, should not exceed 1
    eng. beta_probe = param_input.beta_probe;	% probe step size, larger == faster convergence, smaller == more robust, should not exceed 1
    eng. delta_p = param_input.delta_p;     % LSQ dumping constant, 0 == no preconditioner, 0.1 is usually safe, Preconditioner accelerates convergence and ML methods become approximations of the second order solvers 
    eng. momentum = param_input.momentum;   % add momentum acceleration term to the MLc method, useful if the probe guess is very poor or for acceleration of multilayer solver, but it is quite computationally expensive to be used in conventional ptycho without any refinement. The momentum method works usually well even with the accelerated_gradients option.  eng.momentum = multiplication gain for velocity, eng.momentum == 0 -> no acceleration, eng.momentum == 0.5 is a good value
    eng. accelerated_gradients_start = inf; % iteration number from which the Nesterov gradient acceleration should be applied, this option is supported only for MLc method. It is very computationally cheap way of convergence acceleration. 

    % DM
    eng. pfft_relaxation = 0.05;          % Relaxation in the Fourier domain projection, = 0  for full projection 
    eng. probe_regularization = 0.1;      % Weight factor for the probe update (inertia)

    % ADVANCED OPTIONS                     See for more details: Odstrčil M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
    % position refinement 
    eng. apply_subpix_shift = true;       % apply FFT-based subpixel shift, it is automatically allowed for position refinement
    eng. probe_position_search = param_input.probe_position_search;      % iteration number from which the engine will reconstruct probe positions, from iteration == probe_position_search, assume they have to match geometry model with error less than probe_position_error_max
    eng. probe_geometry_model = {};  % list of free parameters in the geometry model, choose from: {'scale', 'asymmetry', 'rotation', 'shear'
    if param_input.geometry_model_scale
    	eng. probe_geometry_model{end+1} = 'scale';
    end
    if param_input.geometry_model_asymmetry
    	eng. probe_geometry_model{end+1} = 'asymmetry';
    end
    if param_input.geometry_model_rotation
    	eng. probe_geometry_model{end+1} = 'rotation';
    end
    if param_input.geometry_model_shear
    	eng. probe_geometry_model{end+1} = 'shear';
    end

    eng. probe_position_error_max = inf; % maximal expected random position errors, probe prositions are confined in a circle with radius defined by probe_position_error_max and with center defined by original positions scaled by probe_geometry_model
    eng. apply_relaxed_position_constraint = param_input.apply_geometry_model_constraint; % added by YJ. Apply a relaxed constraint to probe positions. default = true. Set to false if there are big jumps in positions.
    eng. update_pos_weight_every = param_input.update_pos_weight_every; % added by YJ. Allow position weight to be updated multiple times. default = inf: only update once.

    % multilayer extension 
    if strcmp(eng.name, 'GPU_MS')
        eng. delta_z = param_input.delta_z * ones(param_input.Nlayers, 1);                     % if not empty, use multilayer ptycho extension , see ML_MS code for example of use, [] == common single layer ptychography , note that delta_z provides only relative propagation distance from the previous layer, ie delta_z can be either positive or negative. If preshift_ML_probe == false, the first layer is defined by position of initial probe plane. It is useful to use eng.momentum for convergence acceleration 
        eng. regularize_layers = param_input.regularize_layers;           % multilayer extension: 0<R<<1 -> apply regularization on the reconstructed object layers, 0 == no regularization, 0.01 == weak regularization that will slowly symmetrize information content between layers 
        eng. preshift_ML_probe = false;       % multilayer extension: if true, assume that the provided probe is reconstructed in center of the sample and the layers are centered around this position 
        eng. layer4pos = [];  % Added by ZC. speficy which layer is used for position correction ; if empty, then default, ceil(Nlayers/2)
        if param_input.layer4pos > 0
            eng. layer4pos = [param_input.layer4pos];  % Added by ZC. speficy which layer is used for position correction ; if empty, then default, ceil(Nlayers/2)
        end
        eng. init_layer_select = [];          % Added by YJ. Select layers in the initial object for pre-processing. If empty (default): use all layers.
        eng. init_layer_preprocess = '';      % Added by YJ. Specify how to pre-process initial layers
                                              % '' or 'all' (default): use all layers (do nothing)
                                              % 'avg': average all layers 
                                              % 'interp': interpolate layers using spline method. Need to specify desired depths in init_layer_interp
        eng. init_layer_interp = [];          % Specify desired depths for interpolation. The depths of initial layers are [1:Nlayer_init]. If empty (default), no interpolation                 
        eng. init_layer_append_mode = param_input.init_layer_append_mode;     % Added by YJ. Specify how to initialize extra layers
                                              % '' or 'vac' (default): add vacuum layers
                                              % 'edge': append 1st or last layers
                                              % 'avg': append averaged layer
        eng. init_layer_scaling_factor = param_input.init_layer_scaling_factor;   % Added by YJ. Scale all layers. Default: 1 (no scaling). Useful when delta_z is changed
        eng. save_images = {'obj_ph_stack','obj_ph_sum','probe','probe_mag','probe_prop_mag'};
    else
        eng. delta_z = [];                     % if not empty, use multilayer ptycho extension , see ML_MS code for example of use, [] == common single layer ptychography , note that delta_z provides only relative propagation distance from the previous layer, ie delta_z can be either positive or negative. If preshift_ML_probe == false, the first layer is defined by position of initial probe plane. It is useful to use eng.momentum for convergence acceleration 
        eng. regularize_layers = 0;            % multilayer extension: 0<R<<1 -> apply regularization on the reconstructed object layers, 0 == no regularization, 0.01 == weak regularization that will slowly symmetrize information content between layers 
        eng. preshift_ML_probe = true;         % multilayer extension: if true, assume that the provided probe is reconstructed in center of the sample and the layers are centered around this position 
        eng. save_images = {'obj_ph','probe_mag','probe'};
    end
    
    % other extensions 
    eng. background = param_input.background;               % average background scattering level, for OMNI values around 0.3 for 100ms, for flOMNI <0.1 per 100ms exposure, see for more details: Odstrcil, M., et al., Optics letters 40.23 (2015): 5574-5577.
    eng. background_width = inf;           % width of the background function in pixels,  inf == flat background, background function is then convolved with the average diffraction pattern in order to account for beam diversion 
    eng. clean_residua = false;            % remove phase residua from reconstruction by iterative unwrapping, it will result in low spatial freq. artefacts -> object can be used as an residua-free initial guess for netx engine
    eng. diff_pattern_blur = param_input.diff_pattern_blur; % account for blurring in diffraction patterns using a gaussian kernel with the given std. Default=0 (no blurring)

    % wavefront & camera geometry refinement     See for more details: Odstrčil M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
    eng. probe_fourier_shift_search = inf; % iteration number from which the engine will: refine farfield position of the beam (ie angle) from iteration == probe_fourier_shift_search
    eng. estimate_NF_distance = inf;       % iteration number from which the engine will: try to estimate the nearfield propagation distance using gradient descent optimization  
    eng. detector_rotation_search = inf;   % iteration number from which the engine will: search for optimal detector rotation, preferably use with option mirror_scan = true , rotation of the detector axis with respect to the sample axis, similar as rotation option in the position refinement geometry model but works also for 0/180deg rotation shared scans 
    eng. detector_scale_search = inf;      % iteration number from which the engine will: refine pixel scale of the detector, can be used to refine propagation distance in ptycho 
    eng. variable_probe = param_input.variable_probe_modes > 0;           % Use SVD to account for variable illumination during a single (coupled) scan, see for more details:  Odstrcil, M. et al. Optics express 24.8 (2016): 8360-8369.
    eng. variable_probe_modes = param_input.variable_probe_modes;         % OPRP settings , number of SVD modes using to describe the probe evolution. 
    eng. variable_probe_smooth = param_input.variable_probe_smooth;        % OPRP settings , enforce of smooth evolution of the OPRP modes -> N is order of polynomial fit used for smoothing, 0 == do not apply any smoothing. Smoothing is useful if only a smooth drift is assumed during the ptycho acquisition 
    eng. variable_intensity = param_input.variable_intensity;       % account to changes in probe intensity

    % extra analysis
    eng. get_fsc_score = false;            % measure evolution of the Fourier ring correlation during convergence 
    eng. mirror_objects = false;           % mirror objects, useful for 0/180deg scan sharing -> geometry refinement for tomography, works only if 2 scans are provided 

    % custom data adjustments, useful for offaxis ptychography
    eng.auto_center_data = false;           % autoestimate the center of mass from data and shift the diffraction patterns so that the average center of mass corresponds to center of mass of the provided probe 
    eng.auto_center_probe = false;          % center the probe position in real space before reconstruction is started 
    eng.manual_center_probe = [param_input.manual_center_probe_y,...
                                param_input.manual_center_probe_x];          % center the probe position in real space before reconstruction is started 

    eng.custom_data_flip = [param_input.dp_custom_fliplr,...
                            param_input.dp_custom_flipud,...
                            param_input.dp_custom_transpose];         % apply custom flip of the data [fliplr, flipud, transpose]  - can be used for quick testing of reconstruction with various flips or for reflection ptychography 
    eng.apply_tilted_plane_correction = ''; % if any(p.sample_rotation_angles([1,2]) ~= 0),  this option will apply tilted plane correction. (a) 'diffraction' apply correction into the data, note that it is valid only for "low NA" illumination  Gardner, D. et al., Optics express 20.17 (2012): 19050-19059. (b) 'propagation' - use tilted plane propagation, (c) '' - will not apply any correction 

    %added by YJ
    if Nfly > 0
        eng.extension = {'fly_scan'};           %arbitrary fly scan
        eng.flyscan_trajectory = 'line';      % Added by YJ. Specify trajectory type for arbitrary-path fly-scan:
                                              %'line' (default): line scan with big jumps. 
                                              %'continuous': contiuous path. 
                                              %'external': load positions from external files
        eng.Nmodes = Nfly;
    end
    eng.plot_results_every = inf;
    eng.save_results_every = param_input.Niter_save_results_every;
    eng.save_init_probe = param_input.save_init_probe; %save initial probe function in the .mat file
    
    resultDir = strcat(param_input.output_dir_base,sprintf(p.scan_string_format,  p.scan_number));
    if ~isempty(p.scan.roi_label); resultDir = [resultDir,'/roi',p.scan.roi_label]; end
    resultDir = fullfile(resultDir,'/',param_input.output_dir_prefix);
    
    eng.extraPrintInfo = strcat('Scan',num2str(p.scan_number(1)));
    output_dir_suffix = param_input.output_dir_suffix;
    if p.detector.upsampling>0
        output_dir_suffix = sprintf([output_dir_suffix,'_dpUpsample%d'], 2^p.detector.upsampling);
    end
    [eng.fout, p.suffix] = generateResultDir(eng, resultDir, output_dir_suffix);
    [p, ~] = core.append_engine(p, eng);    % Adds this engine to the reconstruction process

    %% Run the reconstruction
    tic
    out = core.ptycho_recons(p);
    try
        data_error = out.error_metric.value(end);
    catch
        data_error = inf;
    end
    toc
end
