% Wrapper function for simulating electron ptychography data and evaluating
% its reconstruction quality
% Can be used with BO-GP to find optimal experimental conditions or reconstruction parameters
% Examples: https://doi.org/10.1038/s41598-022-16041-5

function [recon_score] = sim_electron_ptycho_recon(params, varargin)

    parser = inputParser;
    %experimental parameters
    parser.addParameter('field_of_view',  60 , @isnumeric) %scan field of view in angstroms
    parser.addParameter('scan_step_size',  4 , @isnumeric)  %scan step size in angstroms
    parser.addParameter('dose',  1e6 , @isnumeric)  %total electron dose (e/A^2)
    parser.addParameter('N_dp',  128 , @isnumeric)  %size of experimental diffraction pattern in pixels
    parser.addParameter('N_dp_orig',  1024 , @isnumeric)  %size of true diffraction pattern in pixels
    parser.addParameter('max_position_error',  0 , @isnumeric)  %max scan position error in angstroms
    parser.addParameter('voltage',  80 , @isnumeric)  %beam voltage in kV

    % probe parameters 
    parser.addParameter('alpha_max',  0 , @isnumeric) %aperture size in mrad
    parser.addParameter('df', 0 , @isnumeric ) %defocus in nanometers
    parser.addParameter('cs', 0 , @isnumeric ) %third-order spherical aberration in millimeters 
    parser.addParameter('c5', 0 , @isnumeric ) %fifth-order spherical aberration in meters
    parser.addParameter('c7', 0 , @isnumeric ) %seventh-order spherical aberration in meters
    parser.addParameter('f_a2', 0 , @isnumeric ) %twofold astigmatism in nanometers
    parser.addParameter('theta_a2', 0 , @isnumeric ) %twofold azimuthal orientation in radian
    parser.addParameter('f_a3', 0 , @isnumeric ) %threefold astigmatism in micrometers
    parser.addParameter('theta_a3', 0 , @isnumeric ) %threefold azimuthal orientation in radian

    % reconstruction parameters 
    parser.addParameter('eng_name', 'GPU', @ischar)
    parser.addParameter('Niter', 1000, @isnumeric)
    parser.addParameter('Niter_save_results_every', 1000, @isnumeric)
    parser.addParameter('Nprobe', 1, @isnumeric)
    parser.addParameter('grouping', 100, @isnumeric)
    parser.addParameter('method',  'MLs', @ischar)
    parser.addParameter('momentum',  0, @isnumeric)
    parser.addParameter('opt_errmetric',  'L1', @ischar)
    parser.addParameter('apply_multimodal_update', 0, @islogical)
    parser.addParameter('probe_position_search', inf, @isnumeric)

    parser.addParameter('GPU_list', 1 , @isnumeric ) %GPU IDs for reconstruction

    % evaluation parameters 
    parser.addParameter('metric', 'ssim', @ischar)
    parser.addParameter('crop_x', 0, @isnumeric)
    parser.addParameter('crop_y', 0, @isnumeric)
    parser.addParameter('ground_truth_recon', '', @ischar)

    parser.parse(varargin{:})
    r = parser.Results;

    % load all to the param structure 
    par = params; 
    for name = fieldnames(r)'
        if ~isfield(par, name{1}) || ~ismember(name, parser.UsingDefaults) % prefer values in param structure 
            par.(name{1}) = r.(name{1});
        end
    end
    
    %% check inputs
    assert(~isempty(par.ground_truth_recon), 'Ground truth recon file cannot be empty!')

    %% initialize parameters
    base_path = par.base_path;
    par_sim = {}; % for simulation function
    par_sim.field_of_view = par.field_of_view; %scan field of view in angstroms
    par_sim.scan_step_size = par.scan_step_size; %scan step size in angstroms
    par_sim.dose = par.dose; %total electron dose (e/A^2)
    par_sim.N_dp = par.N_dp; %size of experimental diffraction pattern in pixels
    par_sim.N_dp_orig = par.N_dp_orig; %size of true diffraction pattern in pixels
    par_sim.max_position_error = par.max_position_error; %max scan position error in angstroms
    par_sim.voltage = par.voltage; %beam voltage in kV

    par_sim.alpha_max = par.alpha_max; %mrad
    par_sim.probe_df = par.df*10; %angstrom
    par_sim.probe_c3 = par.cs*1e-3/1e-10; %angstrom
    par_sim.probe_c5 = par.c5/1e-10; %angstrom
    par_sim.probe_c7 = par.c7/1e-10; %angstrom
    par_sim.probe_f_a2 = par.f_a2*10; %angstrom
    par_sim.probe_theta_a2 = par.theta_a2; %rad
    par_sim.probe_f_a3 = par.f_a3*1e-6/1e-10; %angstrom
    par_sim.probe_theta_a3 = par.theta_a3; %rad

    data_path = generate_data_path(par.base_data_path, varargin);
    par_sim.output_path = fullfile(base_path, data_path, 'data1');
    par_sim.dx = par.dx;
    par_sim.object = par.object_true;

    if length(par.GPU_list)>1 %assume parallel processing
        t = getCurrentTask;
        t = t.ID;
        gpu_id = par.GPU_list(t);
    else
        gpu_id = par.GPU_list;
    end
    par_sim.gpu_id = gpu_id;

    %% start simulation
    disp('Simulate diffraction patterns...')

    [~, ~, p] = sim_electron_cbed(par_sim);

    disp('Simulate diffraction patterns...done')

    %% run ptycho recon
    disp('Reconstruction...')
    par_recon = {};
    par_recon.Niter = par.Niter;
    par_recon.Niter_save_results_every =  par.Niter_save_results_every;
    par_recon.Nprobe = par.Nprobe;
    par_recon.grouping = par.grouping;
    par_recon.gpu_id = gpu_id;

    par_recon.method = par.method;
    par_recon.verbose_level = 0;

    par_recon.scan_number = 1;
    par_recon.beam_source = 'electron';

    par_recon.Ndp = par_sim.N_dp;
    par_recon.base_path = fullfile(base_path, data_path,'/');
    par_recon.cen_dp_y = floor(par_recon.Ndp/2)+1;
    par_recon.cen_dp_x = floor(par_recon.Ndp/2)+1;

    par_recon.initial_probe_file = fullfile(par_sim.output_path, 'init_probe.mat');
    par_recon.energy = par_sim.voltage;
    par_recon.dk = p.dk; %calculated in sim_electron_cbed
    par_recon.scan_format = 'data%d';
    par_recon.scan_string_format = 'data%d';
    par_recon.roi_label = '0';

    par_recon.detector_name = 'empad';
    par_recon.data_preparator = 'matlab_aps';

    par_recon.src_positions =  'hdf5_pos';
    par_recon.scan_type = 'default';

    par_recon.output_dir_base = fullfile(base_path, data_path, '/');

    [~, eng, ~] = ptycho_recon(par_recon);
    disp('Reconstruction...done')

    %% evaluate ptycho recon
    disp('Evaluate reconstruction...')
    par_eval = {};
    par_eval.file1 = fullfile(eng.fout, sprintf('Niter%d.mat', par_recon.Niter));
	par_eval.file2 = par.ground_truth_recon;
    
    par_eval.crop_y = par.crop_y;
    par_eval.crop_x = par.crop_x;
    
    par_eval.electron = true;
    par_eval.verbose_level = 0;
    par_eval.metric = par.metric;
    recon_score = compare_two_ptycho_recons(par_eval);

    switch par_eval.metric
        case 'frc_1bit'
        otherwise
            recon_score = 1 - recon_score;
    end
    disp('Evaluate reconstruction...done')
    
end

function [base_path] = generate_data_path(base_path, param_var)

    for i=1:2:length(param_var)
        switch param_var{i}
            case 'alpha_max'
                par_format = '_alpha%0.2fmrad';
            case 'df'
                par_format = '_df%0.2fnm';
            case 'cs'
                par_format = '_cs%0.2fmm';
            case 'c5'
                par_format = '_c5%0.2fm';
            case 'c7'
                par_format = '_c7%0.2fm';
            case 'fa2'
                par_format = '_fa2_%0.2fnm';
            case 'theta_a2'
                par_format = '_theta_a2_%0.2frad';
            case 'fa3'
                par_format = '_fa3_%0.2fum';
            case 'theta_a3'
                par_format = '_theta_a3_%0.2frad';
            case 'Ndp'
                par_format = '_Ndp%d';
            case 'N_dp_orig'
                par_format = '_Ndp_orig%d';
            case 'scan_step_size'
                par_format = '_stepsize%0.1fA';
            case 'dose'
                par_format = '_dose%d';
            otherwise
                par_format = ['_', param_var{i}, '%d'];
        end
        base_path = sprintf([base_path, par_format], param_var{i+1});
    end    
end

