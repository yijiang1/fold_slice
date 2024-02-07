clc
close all
clear

%%
par = {}; %basic parameters
par.base_path = '/your/project/directory/bilayer_MoS2/';
par.field_of_view = 36; %scan field of view in angstroms
par.scan_step_size = 4; %scan step size in angstroms
par.N_dp = 128; %size of the experimental diffraction pattern in pixels
par.N_dp_orig = 1024; %size of the true diffraction pattern in pixels
par.N_dp_crop = 1024; %size of the cropped diffraction pattern in pixels
par.max_position_error = 0; %max scan position errors in angstroms
par.voltage = 80;

par.Niter = 500;
par.Niter_save_results_every = par.Niter;
par.grouping = inf;
par.method = 'MLc';
par.Nprobe = 1;

% for evaluation
par.ground_truth_recon = fullfile(pwd, 'utils_electron', 'Bilayer_MoS2_30deg_0p125_2048x2048_ptycho_recon.mat');

par.crop_y = 120;
par.crop_x = 120;

par.overwrite_data = false;
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
par.alpha_max = 20;
par.GPU_list = [1];

if par.dose < inf
    par.base_data_path = strcat('dose', num2str(par.dose));
else
    par.base_data_path = strcat('dose_inf');
end

par.base_data_path = sprintf([par.base_data_path, '_Ndp%d_ss%0.2fA_a%0.2fmrad/probe'], par.N_dp, par.scan_step_size, par.alpha_max);

recon_score = sim_electron_ptycho_recon(par, 'df', 30, 'cs', 0);

%% Use BO-GP to find optimal aberrations
clc; close all; set(0, 'DefaultTextInterpreter', 'none')

par.dose = 1e5;
par.alpha_max = 20;
par.GPU_list = [1];

aberrations = 'df';

bo_verbose = 2;
plot_funcs = {@plotObjectiveModel, @plotMinObjective};

df = optimizableVariable('df', [0, 400]); %nm
cs = optimizableVariable('cs', [0, 2],'Transform','none'); %mm
c5 = optimizableVariable('c5', [0, 1]); %m
c7 = optimizableVariable('c7', [0, 2e3], 'Transform','none'); %m
f_a2 = optimizableVariable('f_a2', [0, 150]); %nm
f_a3 = optimizableVariable('f_a3', [0, 10]); %um
theta_a2 = optimizableVariable('theta_a2', [0, pi]); %rad
theta_a3 = optimizableVariable('theta_a3', [0, 2*pi/3]); %rad

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
par.base_data_path = sprintf([par.base_data_path, '_Ndp%d_ss%0.1fA_a%0.1fmrad/probe'], par.N_dp, par.scan_step_size, par.alpha_max);
par.base_data_path = [par.base_data_path, '_', aberrations, '/probe'];

% Prepare simulation function
switch aberrations
    case 'df'
        func = @(x)sim_electron_ptycho_recon(par, 'df', x.df);
        func_inputs = [df];
        N_init = 5;
        N_eval = 40;
    case 'cs'
        func = @(x)sim_electron_ptycho_recon(par, 'cs', x.cs);
        func_inputs = [cs];
        N_init = 5;
        N_eval = 40;
    case 'c5'
        func = @(x)sim_electron_ptycho_recon(par, 'c5', x.c5);
        func_inputs = [c5];
        N_init = 5;
        N_eval = 40;
    case 'c7'
        func = @(x)sim_electron_ptycho_recon(par, 'c7', x.c7);
        func_inputs = [c7];
        N_init = 5;
        N_eval = 40;
    case 'df_cs'
        func = @(x)sim_electron_ptycho_recon(par, 'df', x.df, 'cs', x.cs);
        func_inputs = [df, cs];
        N_init = 20;
        N_eval = 100;
    case 'f_a2_theta_a2'
        func = @(x)sim_electron_ptycho_recon(par, 'f_a2', x.f_a2, 'theta_a2', x.theta_a2);
        func_inputs = [f_a2, theta_a2];
        N_init = 20;
        N_eval = 100;
   case 'f_a3_theta_a3'
        func = @(x)sim_electron_ptycho_recon(par, 'f_a3', x.f_a3, 'theta_a3', x.theta_a3);
        func_inputs = [f_a3, theta_a3];
        N_init = 20;
        N_eval = 100;
    case 'df_cs_c5'
        func = @(x)sim_electron_ptycho_recon(par, 'df', x.df, 'cs', x.cs, 'c5', x.c5);
        func_inputs = [df, cs, c5];
        N_init = 30;
        N_eval = 200;
end

% Begin BO
results = bayesopt(func, func_inputs,...
    'Verbose', bo_verbose,...
    'AcquisitionFunctionName', 'expected-improvement-plus',...
    'IsObjectiveDeterministic', false,...
    'MaxObjectiveEvaluations', N_eval,...
    'NumSeedPoints', N_init,...
    'PlotFcn',plot_funcs,'UseParallel', N_workers>1);

% Save BO result
save_path = sprintf([par.base_path, '/summary/bo_dose%d_Ndp%d_ss%0.1fA_a%0.1fmrad/'], par.dose, par.N_dp, par.scan_step_size, par.alpha_max);
save_name = [aberrations, '_.mat'];

if ~exist(save_path, 'dir'); mkdir(save_path); end
save(fullfile(save_path, save_name), 'results')

delete(gcp('nocreate'))
disp('Optimization complete')

