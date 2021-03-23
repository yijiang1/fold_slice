% INITIALIZE generate list of default parameters 
% [param] = get_defaults 
% 
%
% returns: 
% ++ param       structure containing parameters for the engines 

function [param] = get_defaults

    %%%%%%%%%%%%%% GPU SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%
    param.use_gpu = true;        % use GPU if possible 
    param.keep_on_gpu = true;    % keep the data all the time on GPU
    param.compress_data = true;  % apply online compress on the GPU data 
    param.gpu_id = []; % default GPU id, [] means choosen by matlab
    param.check_gpu_load = true;
    
    %% basic recontruction parameters 
    %% PIE 
    param.beta_object = 1;
    param.beta_probe = 1;  % step size, faster convergence , more instable ?? 
    %% DM
    param.pfft_relaxation = 0.1; 
    param.probe_inertia = 0.3; % add inertia to the probe reconstruction to avoid oscilations 
    %% general 
    param.share_probe = true;
    param.share_object = false;
    param.delta = 0;  % press values to zero out of the probe area !!  illim < max*delta is removed 
    param.relax_noise = 0.0;  % relaxation for noise, lower => slower convergence, more robust 
    param.positivity_constraint_object = 0; % enforce weak positivity in object 
    param.Nmodes = 1;  %  number of multi apertures , always better to start wih one !! 
    param.probe_modes = 1; % number of probes 
    param.object_modes = 1;  %  number of multi apertures , always better to start wih one !! 
    param.probe_change_start = 1;  % iteration when the probe reconstruction is started
    param.object_change_start = 1; % iteration when the object reconstruction is started
    param.number_iterations = 300;
    param.time_limit = inf; % added by YJ: set a time limit (in seconds) for reconstruction
    param.grouping = inf;
    param.method = 'MLs';
    param.likelihood = 'L1' ; % l1 or poisson,   - choose which likelihood should be used for solver, poisson is suported only for PIE 
    param.verbose_level = 1;
    param.plot_results_every = 50;

    param.remove_residues = false; % autodetect and remove phase residua 
    param.extension = '';

    %% data handling 
    param.upsampling_data_factor = 0;           % assume that the data were created by upsampling using function utils.unbinning 

    param.damped_mask = 5e-3;  % if damped_mask = 0 -> do nothing, if 1>x>0  ->  push masked regions weakly towards measured magnitude value in each iteration
    
    param.background_detection = false; 
    param.background_width = inf;

    %% ADVANCED OPTIONS   
    
    param.object_regular =  [0, 0]; %  enforce smoothness !!!, use between [0-0.1 ]
    param.remove_object_ambiguity = true;    % remove intensity ambiguity between the object and the probes 
    param.variable_probe = false;           % Use SVD to account for variable illumination during a single (coupled) scan
    param.apply_subpix_shift = false;       % apply FFT-based subpixel shift, important for good position refinement but it is slow

    param.probe_geometry_model = {'scale', 'asymmetry', 'rotation', 'shear'};  % list of free parameters in the geometry model
    param.probe_position_search = inf;
    param.apply_relaxed_position_constraint = true; %added by YJ: allow position update without geom model constraint
    param.update_pos_weight_every = inf; %added by YJ: allow position weight to be updated multiple times. Default = inf: only calculate once
    param.max_pos_update_shift = 0.1; %added by YJ: allow user to specify the maximum position update allowed in each iteration. Default = 0.1 (pixel).
    param.probe_position_search_momentum = 0; % added by YJ: enable momentum acceleration for position correction. Default = 0: no acceleration.

    param.probe_fourier_shift_search = inf; 
    param.estimate_NF_distance = inf;
    param.detector_rotation_search = inf;   % rotation of the detector axis with respect to the sample axis, similar as rotation option in the position refinement geometry model but works also for 0/180deg rotation shared scans 
    param.detector_scale_search = inf;      % pixel scale of the detector, can be used to refine propagation distance in ptycho 

    param.apply_multimodal_update = false; % use thibault modes to get higher signal, it can cause isses, not real gain  if blur method is used 
    param.probe_backpropagate = 0; 
    param.beta_LSQ = 0.9;       % use predictive step length
    param.delta_p = 0.1;     % LSQ damping constant 
    param.variable_probe_modes = 1; % OPRP settings 
    param.variable_probe_smooth = 0;% OPRP settings 
    param.variable_intensity = false; % account fort variable intensity
    param.relaxed_object_constrain = 0; % enforce known object (inputs.object_orig)
    param.probe_position_error_max = 10e-9; % max expected error of the stages 
    param.probe_fourier_shift_search = inf; 
    param.momentum = 0;             % use mementume accelerated gradient decsent method 

    param.regularize_layers = 0;    % 0<R<1 -> apply regularization on the reconstructed layers 
    param.preshift_ML_probe = true; % multilayer ptycho extension: if true, assume that the provided probe is reconstructed in center of the sample. 
    
    param. initial_probe_rescaling = true;  % find the optimal scaling correction for the provided probe guess in the initial iteration 
    param. accelerated_gradients_start = inf;  % use accelerated gradients to speed up the convergence
    param. align_shared_objects = false;      % align multiple objects from various scans 

    % extra analysis
    param. get_fsc_score = false;         % measure evolution of the Fourier ring correlation during convergence 
    param. mirror_objects = false;        % mirror objects, useful for 0/180deg scan sharing 
    param. align_shared_objects = false;   % align the objects before sharing them onto single one 

    % fly scans 
    param.flyscan_offset = 0; 
    param.flyscan_dutycycle = 1;
    param.flyscan_intensity = 'varying';  % Added by YJ. Specify how to combine flyscan modes: 'varying' (default) or 'constant'.
    param.flyscan_trajectory = 'line';    % Added by YJ. Specify trajectory type for arbitrary-path fly-scan:
                                          %'line' (default): line scan with big jumps. 
                                          %'continuous': contiuous path. 
                                          %'external': load positions from external files
                                          
    % regularizations added by YJ
    % Remove grid artifacts in the phase image via a windowed Fourier filter (assume raster scan)
    % Based on the idea in https://doi.org/10.1063/1.4993744
    param.rm_grid_artifact_step_size = [0,0];       % scan step size in the [horizontal, vertical] directions. No filter if any of them is 0 (default).
    param.rm_grid_artifact_window_size = [5,5];     % window size in the [horizontal, vertical] directions
    param.rm_grid_artifact_direction = 'xy';        % filter directions: 'x' (horizontal), 'y' (vertical), or 'xy' (default)

    rng('default');
    rng('shuffle');
    
    % convergence check - stop reconstruction if fourier error is larger than the previous one by given (relative) threshold. 
    param.fourier_error_threshold = inf; % default: no convergence check.
    
    % I/O
    param.save_init_probe = false;       % Added by YJ. If true, save initial probe function in the .mat output file. Default is false.
	param.save_images = {'obj_ph','probe'}; % Added by YJ. Save intermediate results as tiff images. Options: {'obj_ph','obj_ph_sum','probe_mag','probe'}
end
