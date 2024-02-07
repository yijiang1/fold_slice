clc
close all
clear

%%
par = {}; %basic parameters
par.base_path = '/your/project/directory/bilayer_MoS2/';
par.field_of_view = 36; %scan field of view in angstroms
par.N_dp = 128; %size of the experimental diffraction pattern in pixels
par.N_dp_orig = 1024; %size of the true diffraction pattern in pixels
par.N_dp_crop = 1024; %size of the cropped diffraction pattern in pixels
par.max_position_error = 0; %max scan position errors in angstroms
par.voltage = 80;

par.Niter = 500;
par.Niter_save_results_every = par.Niter;
par.grouping = inf;
par.method = 'MLc';
par.Nprobe = 2;

% for evaluation
par.ground_truth_recon = fullfile(pwd, 'utils_electron', 'Bilayer_MoS2_30deg_0p125_2048x2048_ptycho_recon.mat');

par.crop_y = 120;
par.crop_x = 120;

par.overwrite_data = true;
par.data_path_format =  'exp_params';

%% load true object
disp('Load test object...')
addpath(fullfile(pwd, 'utils_electron'))
load('Bilayer_MoS2_30deg_0p125_2048x2048_alt.mat')
proj_potentials = padarray(sigmaV, [1024, 1024], 'circular', 'both');
%create a complex object
par.object_true = ones(size(proj_potentials)).*exp(1i*proj_potentials);
par.dx = 0.125; %real-space pixel size in angstrom
disp('Load test object...done')

%% Test one simulation
par.dose = 1e5;
par.GPU_list = [1];

if par.dose < inf
    par.base_data_path = strcat('dose', num2str(par.dose));
else
    par.base_data_path = strcat('dose_inf');
end
par.base_data_path = sprintf([par.base_data_path, '_Ndp%d/probe'], par.N_dp);

recon_score = sim_electron_ptycho_recon(par, 'alpha_max', 10, 'df', 40, 'scan_step_size', 4);

%% Use BO-GP to find optimal aberrations
clc; close all; set(0, 'DefaultTextInterpreter', 'none')

par.dose = 1e5;
%par.alpha_max = 20;
par.GPU_list = [1];

bo_verbose = 2;
plot_funcs = {@plotObjectiveModel, @plotMinObjective};

df = optimizableVariable('df', [0, 200]); %nm
alpha_max = optimizableVariable('alpha_max', [5, 30]); %mrad
step_size = optimizableVariable('step_size', [0.5, 6]); %rad

N_workers = length(par.GPU_list);
if N_workers>1
    delete(gcp('nocreate'))
    c = parcluster('local');
    c.NumWorkers = N_workers;
    p = parpool(c);
end

if par.dose < inf
    par.base_data_path = strcat('dose', num2str(par.dose));
else
    par.base_data_path = strcat('dose_inf');
end
par.base_data_path = sprintf([par.base_data_path, '_Ndp%d/probe'], par.N_dp);

% Prepare simulation function
func = @(x)sim_electron_ptycho_recon(par, 'alpha_max', x.alpha_max, 'df', x.df, 'scan_step_size', x.step_size);
func_inputs = [alpha_max, df, step_size];
N_init = 20;
N_eval = 200;

% Begin BO
results = bayesopt(func, func_inputs,...
    'Verbose', bo_verbose,...
    'AcquisitionFunctionName','expected-improvement-plus',...
    'IsObjectiveDeterministic',false,...
    'MaxObjectiveEvaluations', N_eval,...
    'NumSeedPoints', N_init,...
    'PlotFcn',plot_funcs,'UseParallel', N_workers>1);

%% Save BO result
save_path = sprintf([par.base_path, '/summary/bo_dose%d_Ndp%d/'], par.dose, par.N_dp);
save_name = ['bo_alpha_ss_df.mat'];

if ~exist(save_path, 'dir'); mkdir(save_path); end
save(fullfile(save_path, save_name), 'results')

delete(gcp('nocreate'))
disp('Optimization complete')

