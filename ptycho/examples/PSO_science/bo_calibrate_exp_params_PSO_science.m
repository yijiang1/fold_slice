%clear variables
addpath(strcat(pwd,'/utils/'))
addpath(core.find_base_package)
utils.ccc
%% Step 0: Run the prepare_data script to generate data for ptycho reconstruction
%% Step 1: Prepare data and reconstruction parameters
par = {};
par.verbose_level = 2;
par.scan_number = 1;
par.beam_source = 'electron';

par.Ndp = 256;
par.base_path = '/your/project/directory/PrScO3/science/';

par.scan_format = '%01d';
par.scan_string_format = '%01d';
par.roi_label = '0_Ndp256';

par.cen_dp_y = floor(par.Ndp/2)+1;
par.cen_dp_x = floor(par.Ndp/2)+1;

par.energy = 300;
par.rbf = 26;

par.scan_nx = 64;
par.scan_ny = 64;

par.scan_step_size_x = 0.41;
par.scan_step_size_y = 0.41;

par.scan_custom_fliplr = 1;
par.scan_custom_flipud = 1;
par.scan_custom_transpose = 1;

par.detector_name = 'empad';
par.data_preparator = 'matlab_aps';
par.src_positions =  'matlab_pos';
par.scan_type = 'raster';

par.use_model_probe = true;
par.normalize_init_probe = true;

par.output_dir_base = par.base_path;
par.Niter = 20;
par.Niter_save_results_every = par.Niter;

par.eng_name = 'GPU_MS';
par.method = 'MLs';
par.momentum = 0;

par.Nprobe = 8;
par.grouping = 128;
par.apply_multimodal_update = false;

par.Nlayers = 10;
par.regularize_layers = 1;
par.variable_probe_modes = 1;
par.Ndp_presolve = 128;

%% Step 1.5 (optional): Run a single reconstruction to check parameters
par.GPU_list = 1;
par.rot_ang = 0;
par.alpha_max = 21.4;

defocus = -200;
thickness = 210; %in angstroms

par.output_dir_suffix_base = '';
data_error = ptycho_recon_exp_data(par, 'defocus', defocus, 'thickness', thickness);

%% Step 2: Use Bayesian optimization with Gaussian processes to find experimental parameters that minimize data error
% Note: Parallel BO is generally recommended for multislice reconstructions
close all
par.rot_ang = 0;
par.alpha_max = 21.4;
par.GPU_list = [1,2,3,4,5];

defocus = optimizableVariable('defocus', [-1000, 1000]); %angstroms
thickness = optimizableVariable('thickness', [50, 400]); %angstroms

N_workers = length(par.GPU_list);
if N_workers>1
    delete(gcp('nocreate'))
    c = parcluster('local');
    c.NumWorkers = N_workers;
    p = parpool(c);
end

fun = @(x)ptycho_recon_exp_data(par, 'thickness', x.thickness, 'defocus', x.defocus);
results = bayesopt(fun, [thickness, defocus],...
    'Verbose', 2,...
    'AcquisitionFunctionName', 'expected-improvement-plus',...
    'IsObjectiveDeterministic', false,...
    'MaxObjectiveEvaluations', 50,...
    'NumSeedPoints', N_workers,...
    'PlotFcn', {@plotObjectiveModel, @plotMinObjective}, 'UseParallel', N_workers>1);

delete(gcp('nocreate'))
