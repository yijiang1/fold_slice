clear; close all; clc
addpath('utils')
addpath(find_base_package())
utils.verbose(struct('prefix', 'run_tests'))

%% SELECTED TESTS TO BE RUN 
test_astra = true;    % run tests of ASTRA projectiors in +astra and multiGPU astra wrappers in +tomo, check that they all provide consistent results 
test_mex_function = true;  % test MEX functions used by tomo.block_fun, tomo.Ax_sup_partial, tomo.Axt_sup_partial
test_simulated_data = true; % test uses synthetic data to run all templates under well controlled conditions 
test_real_data = true; % test uses real data to run all templates under real experimental conditions 

%% SELECT BASIC PARAMETERS 
GPU_id = 1; 
base_path = './';           % path to store temporal data


gpuDevice(GPU_id);
utils.report_GPU_usage
pause(1)

%%  TEST ASTRA WRAPPER
% This script will run some basic features in the +astra/ ASTRA wrapper code
% ie reconstruciton , splitting on GPU, splitting between multiple workers 
% compare results with expected values to detect inconsistencies 
if test_astra
    utils.verbose(struct('prefix', 'Test ASTRA'))
    run('tests/astra_wrappers_tests.m')
end
%%  TEST MEX FUNCTION  
% the MEX functions in +utils/private use OpenMP to accelerate memory transfer from 
% large array into small sub array and back
if test_mex_function
    utils.verbose(struct('prefix', 'Test MEX'))
    run('tests/test_MEX_functions.m')
end
clearvars -except GPU_id base_path test_real_data test_simulated_data


%% TEMPLATE FOR AUTOMATIC TOMOGRAPHY CODE TESTS 
% test uses synthetic data to run all templates under well controlled
% conditions, use setting if you want to add noise or other difficulties to
% be tested during alignment 
if test_simulated_data
    utils.verbose(struct('prefix', 'Test simulated data'))
    run('tests/test_tomo_simulated_data.m')
end
clearvars -except GPU_id base_path test_real_data 

if test_real_data
    utils.verbose(struct('prefix', 'Test real data'))
    run('tests/test_tomo_real_data.m')
end

