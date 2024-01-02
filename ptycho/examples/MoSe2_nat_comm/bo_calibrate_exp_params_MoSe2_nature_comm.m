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

par.Ndp = 128;
par.base_path = '/your/project/directory/mose2/nat_comm/';

par.scan_format = 'scan%d';
par.scan_string_format = 'scan%d';
par.roi_label = '0_Ndp128';

par.cen_dp_y = floor(par.Ndp/2)+1;
par.cen_dp_x = floor(par.Ndp/2)+1;

par.energy = 80;
par.dk = 0.0197;

par.scan_nx = 60;
par.scan_ny = 60;

par.scan_step_size_x = 0.85;
par.scan_step_size_y = 0.85;

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
par.Niter_save_results_every =  par.Niter;

par.eng_name = 'GPU';
par.method = 'MLs';
par.momentum = 0;

par.Nprobe = 2;
par.grouping = 90;
par.apply_multimodal_update = true;

%% Step 1.5 (optional): Run a single reconstruction to check parameters
clc
par.GPU_list = 2;
defocus = -450;
par.rot_ang = 30;
par.alpha_max = 21.4;
par.thickness = 0; %single-slice ptycho
par.output_dir_suffix_base = sprintf('_rot_ang%0.1f', par.rot_ang);
data_error = ptycho_recon_exp_data(par, 'defocus', defocus);

%% Step 2.1: Use Bayesian optimization with Gaussian processes to find experimental parameters that minimize data error 
close all
clc
par.rot_ang = 30;
par.alpha_max = 21.4;
defocus = optimizableVariable('defocus', [-1000, 1000]); %angstroms
par.GPU_list = [1];
par.thickness = 0; %single-slice ptycho
par.output_dir_suffix_base = sprintf('_rot_ang%0.1f', par.rot_ang);

N_workers = length(par.GPU_list);
if N_workers>1
    delete(gcp('nocreate'))
    c = parcluster('local');
    c.NumWorkers = N_workers;
    p = parpool(c);
end

fun = @(x)ptycho_recon_exp_data(par, 'defocus', x.defocus);
results = bayesopt(fun, [defocus],...
    'Verbose', 2,...
    'AcquisitionFunctionName', 'expected-improvement-plus',...
    'IsObjectiveDeterministic', false,...
    'MaxObjectiveEvaluations', 20,...
    'NumSeedPoints', 3,...
    'PlotFcn', {@plotObjectiveModel, @plotMinObjective}, 'UseParallel', N_workers>1);

delete(gcp('nocreate'))

%% Step 2.2: Use Bayesian optimization with Gaussian processes to find experimental parameters that minimize data error 
close all
clc
par.alpha_max = 21.4;
par.GPU_list = [1,2,3,4,5];
par.thickness = 0;
par.output_dir_suffix_base = '';

rot_ang = optimizableVariable('rot_ang', [-90, 90]); %degrees 
defocus = optimizableVariable('defocus', [-1000, 1000]); %angstroms

N_workers = length(par.GPU_list);
if N_workers>1
    delete(gcp('nocreate'))
    c = parcluster('local');
    c.NumWorkers = N_workers;
    p = parpool(c);
end

fun = @(x)ptycho_recon_exp_data(par, 'defocus', x.defocus, 'rot_ang', x.rot_ang);
results = bayesopt(fun, [defocus, rot_ang],...
    'Verbose', 2,...
    'AcquisitionFunctionName', 'expected-improvement-plus',...
    'IsObjectiveDeterministic', false,...
    'MaxObjectiveEvaluations', 50,...
    'NumSeedPoints', 5,...
    'PlotFcn', {@plotObjectiveModel, @plotMinObjective}, 'UseParallel', N_workers>1);

delete(gcp('nocreate'))

