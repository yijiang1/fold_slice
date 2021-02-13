% initialize reconstruction parameters for simulated ptycho data
addpath(fullfile(pwd,'/utils/'))
addpath(core.find_base_package)
testmode = false;       % In test mode, lock files are not written and plots are generated every iteration

%% %%%%%%%%%%%%%%%%%% load data parameters %%%%%%%%%%%%%%%%%%%%
initial_probe_file = strcat(base_path_ptycho,'/init_probe.mat');
load(initial_probe_file)

energy = p.voltage;
dk = p.dk;
scan_string_format = 'data%d';
roi_label = '0';

%% %%%%%%%%%%%%%%%%%% recon parameters %%%%%%%%%%%%%%%%%%%%
gpu_id = 1;
Niter = 200;
Nprobe = 3;
var_probe_modes = 0;
%grouping = 20;
Ndp = 128;        %size of diffraction patterns
cen_dp = floor(Ndp/2)+1;
Niter_refined = 0; % more refined reconstruction
get_fsc_score = true;
Niter_save_results = 100;
Niter_plot_results = 20; % set to inf if you don't want any plots
%% General
p = struct();
p.   verbose_level = 2;                            % verbosity for standard output (0-1 for loops, 2-3 for testing and adjustments, >= 4 for debugging)
p.   use_display = Niter_plot_results<Niter;                                      % global switch for display, if [] then true for verbose > 1
p.   scan_number = [1,2];                                    % Multiple scan numbers for shared scans

% Geometry
p.   z = 1;                                             % Distance from object to detector 
p.   asize = [Ndp, Ndp];                                     % Diffr. patt. array size
p.   ctr = [cen_dp, cen_dp];                                       % Diffr. patt. center coordinates (y,x) (empty means middle of the array); e.g. [100 207;100+20 207+10];
p.   beam_source = 'electron';                         % Added by YJ for electron pty. Use relativistic corrected formula for wavelength
p.   dk = dk;                                          % Added by YJ. For determinting pixel size in electron pty.
p.   prop_regime = 'farfield';                              % propagation regime: nearfield, farfield (default), !! nearfield is supported only by GPU engines 
p.   focus_to_sample_distance = [];                         % sample to focus distance, parameter to be set for nearfield ptychography, otherwise it is ignored 
p.   FP_focal_distance = [];                                %  if nonempty -> assume Fourier ptychography configuration, FP_focal_distance = focal length of objective lens for Fourier Ptychography only,
p.   angular_correction_setup = 'none';                         % if src_positions=='orchestra', choose angular correction for specific cSAXS experiment: 'flomni', 'omny', 'lamni', 'none', 
p.   energy = energy;                                           % Energy (in keV), leave empty to use spec entry mokev
p.   sample_rotation_angles = [0,0,0];                      % Offaxis ptychography correction , 3x1 vector rotation around [X,Y,beam] axes in degrees , apply a correction accounting for tilted plane oR the sample and ewald sphere curvature (high NA correction)

%p.   affine_angle = 0;                                     % Not used by ptycho_recons at all. This allows you to define a variable for the affine matrix below and keep it in p for future record. This is used later by the affine_matrix_search.m script
p.   affine_matrix = [1 , 0; 0, 1] ; % Applies affine transformation (e.g. rotation, stretching) to the positions (ignore by = []). Convention [yn;xn] = M*[y;x]. For flOMNI we found in June 2019: = [1 , 0.0003583 ; 5.811e-05 , 1 ]; for OMNY we found in October 2018: = [1 0;tan(0.4*pi/180) 1]; laMNI in June 2018  [1,0.0154;-0.0017,1.01]; laMNI in August [1.01 0.0031; -0.0018 1.00] 

% Scan meta data
p.   src_metadata = 'none';                                 % source of the meta data, following options are supported: 'spec', 'none' , 'artificial' - or add new to +scan/+meta/

% Scan queue
p.   queue.name = '';                                       % specify file queue; currently only 'filelist' is supported
p.   queue.path=[''];      % Folder where the queue of files is defined, note the content of files can overwrite some parameters in p-structure

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
p.   detector.name = 'empad';                           % see +detectors/ folder 
p.   detector.check_2_detpos = [];                          % = []; (ignores)   = 270; compares to dettrx to see if p.ctr should be reversed (for OMNY shared scans 1221122), make equal to the middle point of dettrx between the 2 detector positions
p.   detector.data_prefix = '';                             % Default using current eaccount e.g. e14169_1_
p.   detector.binning = false;                              % = true to perform 2x2 binning of detector pixels, for binning = N do 2^Nx2^N binning
p.   detector.upsampling = false;                           % upsample the measured data by 2^data_upsampling, (transposed operator to the binning), it can be used for superresolution in nearfield ptychography or to account for undersampling in a far-field dataset
p.   detector.burst_frames = 1;                             % number of frames collected per scan position

p.   prepare.data_preparator = 'matlab_aps';                    % data preparator; 'python' or 'matlab' or 'matlab_aps'
p.   prepare.auto_prepare_data = true;                      % if true: prepare dataset from raw measurements if the prepared data does not exist
p.   prepare.force_preparation_data = true;                 % Prepare dataset even if it exists, it will overwrite the file % Default: @prepare_data_2d
p.   prepare.store_prepared_data = false;                    % store the loaded data to h5 even for non-external engines (i.e. other than c_solver)
p.   prepare.prepare_data_function = '';                    % (used only if data should be prepared) custom data preparation function handle;
p.   prepare.auto_center_data = false;                      % if matlab data preparator is used, try to automatically center the diffraction pattern to keep center of mass in center of diffraction

% Scan positions
p.   src_positions = 'hdf5_pos';                           % 'spec', 'orchestra', 'load_from_file', 'matlab_pos' (scan params are defined below) or add new position loaders to +scan/+positions/
p.   positions_file = [''];    %Filename pattern for position files, Example: ['../../specES1/scan_positions/scan_%05d.dat']; (the scan number will be automatically filled in)

p.   spec.motor.fine_motors = {};                           % Y and X motor name for positions, leave empty for defaults
p.   spec.motor.fine_motors_scale = [];                     % ptycho expects real positions in m; 
p.   spec.motor.coarse_motors = {};                         % Coarse sample position for shared object, use {X-motor, Y-motor} 
p.   spec.motor.coarse_motors_scale = [];                   % Scale of the coarse motors (to scale the provided values to meters)

% scan parameters for option src_positions = 'matlab_pos' or 'hdf5_pos';
p.   scan.type = 'default';                                  % {'round', 'raster', 'round_roi', 'custom'}
p.   scan.roi_label = roi_label;                            % For APS data
p.   scan.format = scan_string_format;                      % For APS data format for scan directory generation
%%% PSI %%% 
p.   scan.radius_in = 0;                                    % round scan: interior radius of the round scan
p.   scan.radius_out = 5e-6;                                % round scan: exterior radius of the round scan
p.   scan.nr = 10;                                          % round scan: number of intervals (# of shells - 1)
p.   scan.nth = 3;                                          % round scan: number of points in the first shell
p.   scan.lx = 20e-6;                                       % round_roi scan: width of the roi
p.   scan.ly = 20e-6;                                       % round_roi scan: height of the roi
p.   scan.dr = 1.5e-6;                                      % round_roi scan: shell step size
p.   scan.nx = 10;                                          % raster scan: number of steps in x
p.   scan.ny = 10;                                          % raster scan: number of steps in y
p.   scan.step_size_x = 1e-6;                               % raster scan: step size (grid spacing)
p.   scan.step_size_y = 1e-6;                               % raster scan: step size (grid spacing)
p.   scan.step_randn_offset = 0;                            % raster scan: relative random offset from the ideal periodic grid to avoid the raster grid pathology 
p.   scan.b = 0;                                            % fermat: angular offset
p.   scan.n_max = 1e4;                                      % fermat: maximal number of points generated 
p.   scan.step = 0.5e-6;                                      % fermat: step size 
p.   scan.cenxy = [0,0];                                    % fermat: position of center offset 
p.   scan.roi = [];                                         % Region of interest in the object [xmin xmax ymin ymax] in meters. Points outside this region are not used for reconstruction.
                                                            %  (relative to upper corner for raster scans and to center for round scans)    
p.   scan.custom_positions_source = '';                     % custom: a string name of a function that defines the positions; also accepts mat file with entry 'pos', see +scans/+positions/+mat_pos.m
p.   scan.custom_params = [];                               % custom: the parameters to feed to the custom position function.

% I/O1
p.   prefix = '';                                              % For automatic output filenames. If empty: scan number
p.   suffix = '';              % Optional suffix for reconstruction 
p.   scan_string_format = scan_string_format;                  % format for scan string generation, it is used e.g for plotting and data saving 

%%%p.   base_path = '../../';                                  % base path : used for automatic generation of other paths 
p.   base_path = base_path_ptycho;     % base path : used for automatic generation of other paths 

p.   specfile = '';                                         % Name of spec file to get motor positions and check end of scan, defaut is p.spec_file == p.base_path;
p.   ptycho_matlab_path = '';                               % cSAXS ptycho package path
p.   cSAXS_matlab_path = '';                                % cSAXS base package path
p.   raw_data_path{1} = '';                                 % Default using compile_x12sa_filename, used only if data should be prepared automatically
p.   prepare_data_path = '';                                % Default: base_path + 'analysis'. Other example: '/afs/psi.ch/project/CDI/cSAXS_project/analysis2/'; also supports %u to insert the scan number at a later point (e.g. '/afs/psi.ch/project/CDI/cSAXS_project/analysis2/S%.5u')
p.   prepare_data_filename = [];                            % Leave empty for default file name generation, otherwise use [sprintf('S%05d_data_%03dx%03d',p.scan_number(1), p.asize(1), p.asize(2)) p.prep_data_suffix '.h5'] as default 
p.   save_path{1} = '';                                     % Default: base_path + 'analysis'. Other example: '/afs/psi.ch/project/CDI/cSAXS_project/analysis2/'; also supports %u to insert the scan number at a later point (e.g. '/afs/psi.ch/project/CDI/cSAXS_project/analysis2/S%.5u')
p.   io.default_mask_file = '';                             % load detector mask defined in this file instead of the mask in the detector packages, (used only if data should be prepared) 
p.   io.default_mask_type = 'binary';                       % (used only if data should be prepared) ['binary', 'indices']. Default: 'binary' 
p.   io.file_compression = 0;                               % reconstruction file compression for HDF5 files; 0 for no compression
p.   io.data_compression = 3;                               % prepared data file compression for HDF5 files; 0 for no compression
p.   io.load_prep_pos = false;                              % load positions from prepared data file and ignore positions provided by metadata

p.   io.data_descriptor = 'simulation';                     %added by YJ. A short string that describe data when sending notifications 
p.   io.phone_number = ['2178192766'];                      % phone number for sending messages
p.   io.send_failed_scans_SMS = false;                       % send message if p.queue_max_attempts is exceeded
p.   io.send_finished_recon_SMS = false;                    % send message after the reconstruction is completed
p.   io.send_crashed_recon_SMS = false;                     % send message if the reconstruction crashes
p.   io.SMS_sleep = 1800;                                   % max 1 message per SMS_sleep seconds
p.   io.script_name = mfilename;                             % added by YJ. store matlab script name

p.   artificial_data_file = 'template_artificial_data';     % artificial data parameters, set p.src_metadata = 'artificial' to use this template

%% Reconstruction

% Initial iterate object
p.   model_object = true;                                   % Use model object, if false load it from file 
p.   model.object_type = 'rand';                            % specify how the object shall be created; use 'rand' for a random initial guess; use 'amplitude' for an initial guess based on the prepared data

%p.   initial_iterate_object_file{1} = '/mnt/micdata2/velociprobe/2019-2/RAVEN_Pillar/results/Pillar_ML_recon/analysis/S00000-00999/S00239/S00239_128x128_b0__recons.mat';                   %  use this mat-file as initial guess of object, it is possible to use wild characters and pattern filling, example: '../analysis/S%05i/wrap_*_1024x1024_1_recons*'
%p.   initial_iterate_object_file{1} = '/home/beams0/YJIANG/research/algorithm/ptychography/simulation/sim100/scalingFactor1_stepSizeFactor1/MLs_poisson_p2_g40_pc10_vp1_Niter20/Niter20.mat';                   %  use this mat-file as initial guess of object, it is possible to use wild characters and pattern filling, example: '../analysis/S%05i/wrap_*_1024x1024_1_recons*'
p.   initial_iterate_object_file{1} = '';                   %  use this mat-file as initial guess of object, it is possible to use wild characters and pattern filling, example: '../analysis/S%05i/wrap_*_1024x1024_1_recons*'

% Initial iterate probe
p.   model_probe = false;                                    % Use model probe, if false load it from file 
p.   model.probe_is_focused = true;                         % Model probe is focused (false: just a pinhole)
p.   model.probe_central_stop = true;                       % Model central stop
p.   model.probe_diameter = 170e-6;                         % Model probe pupil diameter
p.   model.probe_central_stop_diameter = 60e-6;             % Model central stop diameter
p.   model.probe_zone_plate_diameter = 180e-6;              % Model probe zone plate diameter
p.   model.probe_outer_zone_width = [];                     % Model probe zone plate outermost zone width (not used if not a focused probe) 
p.   model.probe_propagation_dist = 3e-3;                 % Model probe propagation distance (pinhole <-> sample for unfocused, focal-plane <-> sample for focused)
p.   model.probe_focal_length = 51e-3;                      % Model probe focal length (used only if model_is_focused is true
                                                            %   AND model_outer_zone_width is empty)
p.   model.probe_upsample = 10;                             % Model probe upsample factor (for focused probes)

%Use probe from this mat-file (not used if model_probe is true)
p.   initial_probe_file = initial_probe_file; 

p.   probe_file_propagation = 0.0e-3;                            % Distance for propagating the probe from file in meters, = 0 to ignore

% Shared scans - Currently working only for sharing probe and object
p.   share_probe  = 0;                                      % Share probe between scans. Can be either a number/boolean or a list of numbers, specifying the probe index; e.g. [1 2 2] to share the probes between the second and third scan. 
p.   share_object = 0;                                      % Share object between scans. Can be either a number/boolean or a list of numbers, specifying the object index; e.g. [1 2 2] to share the objects between the second and third scan. 

% Modes
p.   probe_modes  = Nprobe;                                 % Number of coherent modes for probe
p.   object_modes = 1;                                      % Number of coherent modes for object
% Mode starting guess
p.   mode_start_pow = [0.02];                               % Normalized intensity on probe modes > 1. Can be a number (all higher modes equal) or a vector
p.   mode_start = 'herm';                                   % (for probe) = 'rand', = 'herm' (Hermitian-like base), = 'hermver' (vertical modes only), = 'hermhor' (horizontal modes only)
p.   ortho_probes = true;                                   % orthogonalize probes after each engine

%% Plot, save and analyze
p.   plot.prepared_data = false;                         % plot prepared data
p.   plot.interval = [];                                    % plot each interval-th iteration, does not work for c_solver code
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
p.   plot.residua = false;                                  % highlight phase-residua in the image of the reconstructed phase

p.   save.external = ~testmode;                             % Use a new Matlab session to run save final figures (saves ~6s per reconstruction). Please be aware that this might lead to an accumulation of Matlab sessions if your single reconstruction is very fast.
p.   save.store_images = false;                              % Write preview images containing the final reconstructions in [p.base_path,'analysis/online/ptycho/'] if p.use_display = 0 then the figures are opened invisible in order to create the nice layout. It writes images in analysis/online/ptycho
p.   save.store_images_intermediate = false;                % save images to disk after each engine
p.   save.store_images_ids = 1:4;                           % identifiers  of the figure to be stored, 1=obj. amplitude, 2=obj. phase, 3=probes, 4=errors, 5=probes spectrum, 6=object spectrum
p.   save.store_images_format = 'png';                      % data type of the stored images jpg or png 
p.   save.store_images_dpi = 150;                           % DPI of the stored bitmap images 
p.   save.exclude = {'fmag', 'fmask', 'illum_sum'};         % exclude variables to reduce the file size on disk
p.   save.save_reconstructions_intermediate = false;        % save final object and probes after each engine
p.   save.save_reconstructions = false;                      % save reconstructed object and probe when full reconstruction is finished 
p.   save.output_file = 'mat';                               % data type of reconstruction file; 'h5' or 'mat'

%% ENGINES
% --------- GPU engines  -------------   See for more details: Odstrčil M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
eng = struct();                        % reset settings for this engine 
eng. name = 'GPU';    
eng. use_gpu = true;                   % if false, run CPU code, but it will get very slow 
eng. keep_on_gpu = true;               % keep data + projections on GPU, false is useful for large data if DM is used
eng. compress_data = false;             % use automatic online memory compression to limit need of GPU memory
eng. gpu_id = gpu_id;                      % default GPU id, [] means choosen by matlab
eng. check_gpu_load = true;            % check available GPU memory before starting GPU engines 

% general
eng. number_iterations = Niter;          % number of iterations for selected method 
eng. asize_presolve = [];      % crop data to "asize_presolve" size to get low resolution estimate that can be used in the next engine as a good initial guess 
%eng. share_probe = 1;                 % Share probe between scans. Can be either a number/boolean or a list of numbers, specifying the probe index; e.g. [1 2 2] to share the probes between the second and third scan. 
%eng. share_object = 0;                % Share object between scans. Can be either a number/boolean or a list of numbers, specifying the object index; e.g. [1 2 2] to share the objects between the second and third scan. 
eng. align_shared_objects = false;     % before merging multiple unshared objects into one shared, the object will be aligned and the probes shifted by the same distance -> use for alignement and shared reconstruction of drifting scans  

eng. method = 'MLc';                   % choose GPU solver: DM, ePIE, hPIE, MLc, Mls, -- recommended are MLc and MLs
eng. opt_errmetric = 'L1' ;            % optimization likelihood - poisson, L1
eng. grouping = grouping;                    % size of processed blocks, larger blocks need more memory but they use GPU more effeciently, !!! grouping == inf means use as large as possible to fit into memory 
                                       % * for hPIE, ePIE, MLs methods smaller blocks lead to faster convergence, 
                                       % * for MLc the convergence is similar 
                                       % * for DM is has no effect on convergence
eng. probe_modes  = p.probe_modes;                % Number of coherent modes for probe
eng. object_change_start = 1;          % Start updating object at this iteration number
eng. probe_change_start = 1;           % Start updating probe at this iteration number

% regularizations
eng. reg_mu = 0;                       % Regularization (smooting) constant ( reg_mu = 0 for no regularization)
eng. delta = 0;                        % press values to zero out of the illumination area in th object, usually 1e-2 is enough 
eng. positivity_constraint_object = 0; % enforce weak (relaxed) positivity in object, ie O = O*(1-a)+a*|O|, usually a=1e-2 is already enough. Useful in conbination with OPRP or probe_fourier_shift_search  

eng. apply_multimodal_update = true; % apply all incoherent modes to object, it can cause isses if the modes collect some crap 
eng. probe_backpropagate = 0;         % backpropagation distance the probe mask, 0 == apply in the object plane. Useful for pinhole imaging where the support can be applied  at the pinhole plane
eng. probe_support_radius = [];       % Normalized radius of circular support, = 1 for radius touching the window    
eng. probe_support_fft = false;       % assume that there is not illumination intensity out of the central FZP cone and enforce this contraint. Useful for imaging with focusing optics. Helps to remove issues from the gaps between detector modules.

% basic recontruction parameters 
% PIE / ML methods                    % See for more details: Odstrčil M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
eng. beta_object = 1;                 % object step size, larger == faster convergence, smaller == more robust, should not exceed 1
eng. beta_probe = 1;                  % probe step size, larger == faster convergence, smaller == more robust, should not exceed 1
eng. delta_p = 0.1;                   % LSQ dumping constant, 0 == no preconditioner, 0.1 is usually safe, Preconditioner accelerates convergence and ML methods become approximations of the second order solvers 
eng. momentum = 0;                    % add momentum acceleration term to the MLc method, useful if the probe guess is very poor or for acceleration of multilayer solver, but it is quite computationally expensive to be used in conventional ptycho without any refinement. 
                                      % The momentum method works usually well even with the accelerated_gradients option.  eng.momentum = multiplication gain for velocity, eng.momentum == 0 -> no acceleration, eng.momentum == 0.5 is a good value
                                      % momentum is enabled only when par.Niter < par.accelerated_gradients_start;
eng. accelerated_gradients_start = inf; % iteration number from which the Nesterov gradient acceleration should be applied, this option is supported only for MLc method. It is very computationally cheap way of convergence acceleration. 

% DM
eng. pfft_relaxation = 0.05;          % Relaxation in the Fourier domain projection, = 0  for full projection 
eng. probe_regularization = 0.1;      % Weight factor for the probe update (inertia)

% ADVANCED OPTIONS                     See for more details: Odstrčil M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
% position refinement 
eng. apply_subpix_shift = true;       % apply FFT-based subpixel shift, it is automatically allowed for position refinement
eng. probe_position_search = inf;      % iteration number from which the engine will reconstruct probe positions, from iteration == probe_position_search, assume they have to match geometry model with error less than probe_position_error_max
%eng. probe_geometry_model = {'rotation','scale','shear','asymmetry'};  % list of free parameters in the geometry model, choose from: {'scale', 'asymmetry', 'rotation', 'shear'}
eng. probe_geometry_model = {};  % list of free parameters in the geometry model, choose from: {'scale', 'asymmetry', 'rotation', 'shear'}
eng. probe_position_error_max = inf; % maximal expected random position errors, probe prositions are confined in a circle with radius defined by probe_position_error_max and with center defined by original positions scaled by probe_geometry_model
eng.apply_relaxed_position_constraint = false; % added by YJ. Apply a relaxed constraint to probe positions. default = true. Set to false if there are big jumps in positions.
eng.update_pos_weight_every = inf; % added by YJ. Allow position weight to be updated multiple times. default = true: only update once.

% multilayer extension
eng. delta_z = [];                     % if not empty, use multilayer ptycho extension , see ML_MS code for example of use, [] == common single layer ptychography , note that delta_z provides only relative propagation distance from the previous layer, ie delta_z can be either positive or negative. If preshift_ML_probe == false, the first layer is defined by position of initial probe plane. It is useful to use eng.momentum for convergence acceleration 
eng. regularize_layers = 0;            % multilayer extension: 0<R<<1 -> apply regularization on the reconstructed object layers, 0 == no regularization, 0.01 == weak regularization that will slowly symmetrize information content between layers 
eng. preshift_ML_probe = true;         % multilayer extension: if true, assume that the provided probe is reconstructed in center of the sample and the layers are centered around this position 

% other extensions 
eng. background = 0;               % average background scattering level, for OMNI values around 0.3 for 100ms, for flOMNI <0.1 per 100ms exposure, see for more details: Odstrcil, M., et al., Optics letters 40.23 (2015): 5574-5577.

eng. background_width = inf;           % width of the background function in pixels,  inf == flat background, background function is then convolved with the average diffraction pattern in order to account for beam diversion 
eng. clean_residua = false;            % remove phase residua from reconstruction by iterative unwrapping, it will result in low spatial freq. artefacts -> object can be used as an residua-free initial guess for netx engine

% wavefront & camera geometry refinement     See for more details: Odstrčil M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
eng. probe_fourier_shift_search = inf; % iteration number from which the engine will: refine farfield position of the beam (ie angle) from iteration == probe_fourier_shift_search
eng. estimate_NF_distance = inf;       % iteration number from which the engine will: try to estimate the nearfield propagation distance using gradient descent optimization  
eng. detector_rotation_search = inf;   % iteration number from which the engine will: search for optimal detector rotation, preferably use with option mirror_scan = true , rotation of the detector axis with respect to the sample axis, similar as rotation option in the position refinement geometry model but works also for 0/180deg rotation shared scans 
eng. detector_scale_search = inf;      % iteration number from which the engine will: refine pixel scale of the detector, can be used to refine propagation distance in ptycho 
eng. variable_probe = var_probe_modes>0;           % Use SVD to account for variable illumination during a single (coupled) scan, see for more details:  Odstrcil, M. et al. Optics express 24.8 (2016): 8360-8369.
eng. variable_probe_modes = var_probe_modes;         % OPRP settings , number of SVD modes using to describe the probe evolution. 
eng. variable_probe_smooth = 0;        % OPRP settings , enforce of smooth evolution of the OPRP modes -> N is order of polynomial fit used for smoothing, 0 == do not apply any smoothing. Smoothing is useful if only a smooth drift is assumed during the ptycho acquisition 
eng. variable_intensity = true;       % account to changes in probe intensity

% extra analysis
eng. get_fsc_score = get_fsc_score;            % measure evolution of the Fourier ring correlation during convergence 
eng. mirror_objects = false;           % mirror objects, useful for 0/180deg scan sharing -> geometry refinement for tomography, works only if 2 scans are provided 

% custom data adjustments, useful for offaxis ptychography
eng.auto_center_data = false;           % autoestimate the center of mass from data and shift the diffraction patterns so that the average center of mass corresponds to center of mass of the provided probe 
eng.auto_center_probe = false;          % center the probe position in real space before reconstruction is started 
eng.custom_data_flip = [0,0,0];         % apply custom flip of the data [fliplr, flipud, transpose]  - can be used for quick testing of reconstruction with various flips or for reflection ptychography 
eng.apply_tilted_plane_correction = ''; % if any(p.sample_rotation_angles([1,2]) ~= 0),  this option will apply tilted plane correction. (a) 'diffraction' apply correction into the data, note that it is valid only for "low NA" illumination  Gardner, D. et al., Optics express 20.17 (2012): 19050-19059. (b) 'propagation' - use tilted plane propagation, (c) '' - will not apply any correction 

%% added by YJ
eng.plot_results_every = Niter_plot_results;
eng.save_results_every = Niter_save_results;
eng.save_phase_image = true;
eng.save_probe_mag = true;
eng.save_sub_objects = false;
eng.avg_photon_threshold = 0; %Added by YJ. Check averaged photon count per pixel during pre-processing. Stop if smaller than the threshold (default = 0.01);

resultDir = p.base_path;
eng.fout =  generateResultDir(eng, resultDir);
[eng.fout, p.suffix] =  generateResultDir(eng, resultDir);
fout_old = eng.fout;
% add engine
[p, ~] = core.append_engine(p, eng);    % Adds this engine to the reconstruction process

if Niter_refined>0
    %%%%%%%%%%%%%%%%%%%%%% refined recon %%%%%%%%%%%%%%%%%%%%%%
    %eng = struct();                        % reset settings for this engine 
    eng. name = 'GPU';    
    eng. use_gpu = true;                   % if false, run CPU code, but it will get very slow 
    eng. keep_on_gpu = true;               % keep data + projections on GPU, false is useful for large data if DM is used
    eng. compress_data = true;             % use automatic online memory compression to limit need of GPU memory
    eng. gpu_id = gpu_id;                      % default GPU id, [] means choosen by matlab
    eng. check_gpu_load = true;            % check available GPU memory before starting GPU engines 

    % general 
    eng. number_iterations = Niter_refined;          % number of iterations for selected method 
    eng. asize_presolve = [];      % crop data to "asize_presolve" size to get low resolution estimate that can be used in the next engine as a good initial guess 
    eng. share_probe = p.share_probe;                 % Share probe between scans. Can be either a number/boolean or a list of numbers, specifying the probe index; e.g. [1 2 2] to share the probes between the second and third scan. 
    eng. share_object = p.share_object;                % Share object between scans. Can be either a number/boolean or a list of numbers, specifying the object index; e.g. [1 2 2] to share the objects between the second and third scan. 
    eng. align_shared_objects = false;     % before merging multiple unshared objects into one shared, the object will be aligned and the probes shifted by the same distance -> use for alignement and shared reconstruction of drifting scans  

    eng. method = 'MLc';                   % choose GPU solver: DM, ePIE, hPIE, MLc, Mls, -- recommended are MLc and MLs
    eng. opt_errmetric = 'poisson' ;            % optimization likelihood - poisson, L1
    eng. grouping = grouping;                   % size of processed blocks, larger blocks need more memory but they use GPU more effeciently, !!! grouping == inf means use as large as possible to fit into memory 
                                           % * for hPIE, ePIE, MLs methods smaller blocks lead to faster convergence, 
                                           % * for MLc the convergence is similar 
                                           % * for DM is has no effect on convergence
    %eng. probe_modes  = 1;                % Number of coherent modes for probe
    eng. object_change_start = 1;          % Start updating object at this iteration number
    eng. probe_change_start = 1;           % Start updating probe at this iteration number
    eng. probe_modes  = p.probe_modes;                % Number of coherent modes for probe

    % regularizations
    eng. reg_mu = 0;                       % Regularization (smooting) constant ( reg_mu = 0 for no regularization)
    eng. delta = 0;                        % press values to zero out of the illumination area in th object, usually 1e-2 is enough 
    eng. positivity_constraint_object = 0; % enforce weak (relaxed) positivity in object, ie O = O*(1-a)+a*|O|, usually a=1e-2 is already enough. Useful in conbination with OPRP or probe_fourier_shift_search  

    eng. apply_multimodal_update = true; % apply all incoherent modes to object, it can cause isses if the modes collect some crap 
    eng. probe_backpropagate = 0;         % backpropagation distance the probe mask, 0 == apply in the object plane. Useful for pinhole imaging where the support can be applied  at the pinhole plane
    eng. probe_support_radius = [];       % Normalized radius of circular support, = 1 for radius touching the window    
    eng. probe_support_fft = false;       % assume that there is not illumination intensity out of the central FZP cone and enforce this contraint. Useful for imaging with focusing optics. Helps to remove issues from the gaps between detector modules.

    % basic recontruction parameters 
    % PIE / ML methods                    % See for more details: Odstrčil M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
    eng. beta_object = 1;                 % object step size, larger == faster convergence, smaller == more robust, should not exceed 1
    eng. beta_probe = 1;                  % probe step size, larger == faster convergence, smaller == more robust, should not exceed 1
    eng. delta_p = 0.1;                   % LSQ dumping constant, 0 == no preconditioner, 0.1 is usually safe, Preconditioner accelerates convergence and ML methods become approximations of the second order solvers 
    eng. momentum = 0;                    % add momentum acceleration term to the MLc method, useful if the probe guess is very poor or for acceleration of multilayer solver, but it is quite computationally expensive to be used in conventional ptycho without any refinement. The momentum method works usually well even with the accelerated_gradients option.  eng.momentum = multiplication gain for velocity, eng.momentum == 0 -> no acceleration, eng.momentum == 0.5 is a good value
    eng. accelerated_gradients_start = inf; % iteration number from which the Nesterov gradient acceleration should be applied, this option is supported only for MLc method. It is very computationally cheap way of convergence acceleration. 

    % DM
    eng. pfft_relaxation = 0.05;          % Relaxation in the Fourier domain projection, = 0  for full projection 
    eng. probe_regularization = 0.1;      % Weight factor for the probe update (inertia)

    % ADVANCED OPTIONS                     See for more details: Odstrčil M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
    % position refinement 
    eng. apply_subpix_shift = true;       % apply FFT-based subpixel shift, it is automatically allowed for position refinement
    eng. probe_position_search = inf;      % iteration number from which the engine will reconstruct probe positions, from iteration == probe_position_search, assume they have to match geometry model with error less than probe_position_error_max
    %eng. probe_geometry_model = {'scale', 'asymmetry', 'rotation', 'shear'};  % list of free parameters in the geometry model, choose from: {'scale', 'asymmetry', 'rotation', 'shear'}
    eng. probe_geometry_model = {};  % list of free parameters in the geometry model, choose from: {'scale', 'asymmetry', 'rotation', 'shear'}
    eng. probe_position_error_max = inf; % maximal expected random position errors, probe prositions are confined in a circle with radius defined by probe_position_error_max and with center defined by original positions scaled by probe_geometry_model
    eng. apply_relaxed_position_constraint = false; % added by YJ. Apply a relaxed constraint to probe positions. default = true. Set to false if there are big jumps in positions.
    eng. update_pos_weight_every = inf; % added by YJ. Allow position weight to be updated multiple times. default = inf: only update once.

    % multilayer extension 
    eng. delta_z = [];                     % if not empty, use multilayer ptycho extension , see ML_MS code for example of use, [] == common single layer ptychography , note that delta_z provides only relative propagation distance from the previous layer, ie delta_z can be either positive or negative. If preshift_ML_probe == false, the first layer is defined by position of initial probe plane. It is useful to use eng.momentum for convergence acceleration 
    eng. regularize_layers = 0;            % multilayer extension: 0<R<<1 -> apply regularization on the reconstructed object layers, 0 == no regularization, 0.01 == weak regularization that will slowly symmetrize information content between layers 
    eng. preshift_ML_probe = true;         % multilayer extension: if true, assume that the provided probe is reconstructed in center of the sample and the layers are centered around this position 

    % other extensions 
    eng. background = 0;               % average background scattering level, for OMNI values around 0.3 for 100ms, for flOMNI <0.1 per 100ms exposure, see for more details: Odstrcil, M., et al., Optics letters 40.23 (2015): 5574-5577.
    eng. background_width = inf;           % width of the background function in pixels,  inf == flat background, background function is then convolved with the average diffraction pattern in order to account for beam diversion 
    eng. clean_residua = false;            % remove phase residua from reconstruction by iterative unwrapping, it will result in low spatial freq. artefacts -> object can be used as an residua-free initial guess for netx engine

    % wavefront & camera geometry refinement     See for more details: Odstrčil M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
    eng. probe_fourier_shift_search = inf; % iteration number from which the engine will: refine farfield position of the beam (ie angle) from iteration == probe_fourier_shift_search
    eng. estimate_NF_distance = inf;       % iteration number from which the engine will: try to estimate the nearfield propagation distance using gradient descent optimization  
    eng. detector_rotation_search = inf;   % iteration number from which the engine will: search for optimal detector rotation, preferably use with option mirror_scan = true , rotation of the detector axis with respect to the sample axis, similar as rotation option in the position refinement geometry model but works also for 0/180deg rotation shared scans 
    eng. detector_scale_search = inf;      % iteration number from which the engine will: refine pixel scale of the detector, can be used to refine propagation distance in ptycho 
    eng. variable_probe = var_probe_modes>0;           % Use SVD to account for variable illumination during a single (coupled) scan, see for more details:  Odstrcil, M. et al. Optics express 24.8 (2016): 8360-8369.
    eng. variable_probe_modes = var_probe_modes;         % OPRP settings , number of SVD modes using to describe the probe evolution. 
    eng. variable_probe_smooth = 0;        % OPRP settings , enforce of smooth evolution of the OPRP modes -> N is order of polynomial fit used for smoothing, 0 == do not apply any smoothing. Smoothing is useful if only a smooth drift is assumed during the ptycho acquisition 
    eng. variable_intensity = var_probe_modes;       % account to changes in probe intensity

    % extra analysis
    eng. get_fsc_score = get_fsc_score;            % measure evolution of the Fourier ring correlation during convergence 
    eng. mirror_objects = false;           % mirror objects, useful for 0/180deg scan sharing -> geometry refinement for tomography, works only if 2 scans are provided 

    % custom data adjustments, useful for offaxis ptychography
    eng.auto_center_data = false;           % autoestimate the center of mass from data and shift the diffraction patterns so that the average center of mass corresponds to center of mass of the provided probe 
    eng.auto_center_probe = false;          % center the probe position in real space before reconstruction is started 
    eng.custom_data_flip = [0,0,0];         % apply custom flip of the data [fliplr, flipud, transpose]  - can be used for quick testing of reconstruction with various flips or for reflection ptychography 
    eng.apply_tilted_plane_correction = ''; % if any(p.sample_rotation_angles([1,2]) ~= 0),  this option will apply tilted plane correction. (a) 'diffraction' apply correction into the data, note that it is valid only for "low NA" illumination  Gardner, D. et al., Optics express 20.17 (2012): 19050-19059. (b) 'propagation' - use tilted plane propagation, (c) '' - will not apply any correction 

    %added by YJ
    eng.plot_results_every = Niter_plot_results;
    eng.save_results_every = Niter_save_results;
    eng.save_phase_image = true;
    eng.save_probe_mag = true;
    [eng.fout, p.suffix] =  generateResultDir(eng, fout_old);

    fout_old = eng.fout;

    eng.extraPrintInfo = strcat('simulation-',num2str(p.scan_number));

    [p, ~] = core.append_engine(p, eng);    % Adds this engine to the reconstruction process
end
%% Run the reconstruction
caller = dbstack;
if length(caller)==1
    tic
    out = core.ptycho_recons(p);
    toc
end
