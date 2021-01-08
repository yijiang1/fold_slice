%% TEMPLATE FOR AUTOMATIC TOMOGRAPHY CODE TESTS 
% perform tests on measured dataset stored in /das/work/p16/p16812/

cd(fullfile( fileparts(mfilename('fullpath')), '..'))
addpath('tests')
addpath('utils')
addpath('./')
addpath(find_base_package)
clearvars -except par0 tested_templates scratch_path GPU_id base_path test_simulated_data test_real_data base_path

 
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Edit this section %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

datasets = 1:5;  % 1-nature chip, 2-FFC particle, 3-retina, 4-local tomo, 5-lamni chip

verbose_level = -1;    % -1 = keep very quiet the reconstructions 

scratch_path = '/das/work/p16/p16812/';  % path to the cSAXS scratch p-folder where are saved the test datasets 


if ~exist('GPU_id', 'var'); GPU_id = [1]; end
if ~exist('base_path', 'var'); base_path = '../'; end
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%



utils.verbose(verbose_level)
setenv('TMP',[base_path,'/tmp']) % TEMP for matlab scripts 

utils.verbose(-1, '=== Searching for data in %s ====', scratch_path)

for dataset = datasets
   
    par0 = struct(); 
    par0.GPU_list = GPU_id; 

switch dataset
    case 1
        %% Nature chip 2016
        % test: large projections, phase residua/errors from sharp transitions 
        tested_templates ={'template_tomo_recons'} ;
        par0.tomo_id = [];
        par0.scanstomo = [2718:3925];
        par0.fileprefix='online_';                  % string at the beginning of the filename, related to reconstruction name
        par0.filesuffix = '_600x600_wrap_1_c';   % string searched at the end of the data filenames,  No need to add _c or _recons, it will look for it 
        par0.file_extension = 'mat';
        par0.analysis_path = fullfile(scratch_path, 'data/e16622_tomo_nature_chip_2016/analysis/');
        par0.surface_calib_file = [];
        par0.omnyposfile = fullfile(scratch_path, 'data/e16622_tomo_nature_chip_2016/specES1/scan_positions/scan_%05d.dat');		%Filename pattern for Orchestra interferometer position files   
        par0.OMNY_angle_file = fullfile(scratch_path, 'data/e16622_tomo_nature_chip_2016/specES1/dat-files/tomography_scannumbers.txt');  % Filename with angles
        par0.max_residua_limit = inf;        % limit used to determine which projection have failed 

    case 2
        %% johanness FCC catalyst
        % test: large projections, low freq. errors, periodic artefacts, vertically asymmetric sample 
        tested_templates ={'template_tomo_recons'} ;
        par0.tomo_id = [];
        par0.scanstomo = [500:1349];
        par0.fileprefix='offline_';                  % string at the beginning of the filename, related to reconstruction name
        par0.filesuffix = '500x500_run_1_recons';   % string searched at the end of the data filenames,  No need to add _c or _recons, it will look for it 
        par0.file_extension = 'mat';
        par0.analysis_path = fullfile(scratch_path, 'data/e16410_tomo_FCC_particle/analysis/');
        par0.surface_calib_file = [];
        par0.omnyposfile = fullfile(scratch_path, 'data/e16410_tomo_FCC_particle/specES1/scan_positions/scan_%05d.dat');		%Filename pattern for Orchestra interferometer position files   
        par0.OMNY_angle_file = fullfile(scratch_path, 'data/e16410_tomo_FCC_particle/specES1/dat-files/tomography_scannumbers.txt');  % Filename with angles
        % solve in lower resolution 
        par0.downsample_projections = 1; 

    case 3
        %% retina from OMNY
        % test: large projections, low freq. errors, huge phase jumps 
        tested_templates ={'template_tomo_recons'} ;
        par0.tomo_id = [];
        par0.scanstomo = [1925:2384];
        par0.fileprefix='online_wrap_';                  % string at the beginning of the filename, related to reconstruction name
        par0.filesuffix = '_452x452_run_1_c';   % string searched at the end of the data filenames,  No need to add _c or _recons, it will look for it 
        par0.file_extension = 'mat';
        par0.analysis_path = fullfile(scratch_path, 'data/e15634_retina_2015_OMNY/analysis/');
        par0.surface_calib_file = [];
        par0.omnyposfile = fullfile(scratch_path, 'data/e15634_retina_2015_OMNY/specES1/omny_recontruct/scan_%05d.dat');		%Filename pattern for Orchestra interferometer position files   
        par0.OMNY_angle_file = fullfile(scratch_path, 'data/e15634_retina_2015_OMNY/specES1/dat-files/omny_scannumbers.txt');  % Filename with angles
        par0.downsample_projections = 1;     % downsample projections by factor of 2^x, set 0 to do nothing and 1,2,.. for different levels of projection binning 
        par0.auto_alignment = false; 
        par0.get_auto_calibration = false; 
        
    case 4
        %% local tomo dataset !! THIS TEST TAKES ~1 HOUR and requires 200GB of RAM !!
        tested_templates ={'template_tomo_interior'} ;
        par0.tomo_id = []; % [68:74]; % Either scan numbers or tomo_id can be given, but not both, if not provided leave tomo_id=[]
        par0.scanstomo = 1700:7690; 
        par0.lowres_tomo_path =fullfile(scratch_path, 'data/e17312_localtomo_FCC_particle/tomogram_delta_S00089_to_S01525_ram-lak_freqscl_1.00.mat'); 

        % IO loading 
        par0.fileprefix='';                  % string at the beginning of the filename, related to reconstruction name
        par0.filesuffix = '_recons';   % string searched at the end of the data filenames,  No need to add _c or _recons, it will look for it 
        par0.file_extension = 'h5';
        par0.downsample_projections = 0;     % downsample projections by factor of 2^x, set 0 to do nothing and 1,2,.. for different levels of projection binning 
        par0.analysis_path = fullfile(scratch_path, 'data/e17312_localtomo_FCC_particle/analysis/');

        par0.clip_amplitude_quantile = 0.95;  %  clip amplitudes in the loaded projections that are exceeding given quantile, if par0.clip_amplitude_quantile == 1, do nothing 
        par0.max_residua_limit = 100;        % limit used to determine which projection have failed 
        par0.surface_calib_file = [];
        par0.omnyposfile = fullfile(scratch_path, 'data/e17312_localtomo_FCC_particle/specES1/scan_positions/scan_%05d.dat');		%Filename pattern for Orchestra interferometer position files   
        par0.OMNY_angle_file = fullfile(scratch_path, 'data/e17312_localtomo_FCC_particle/specES1/dat-files/tomography_scannumbers.txt');  % Filename with angles

        % Other
        par0.save_memory = false;        % try to limit use of RAM 
            par0.inplace_processing = par0.save_memory; % process object_stack using inplace operations to save memory 
            par0.fp16_precision     = true; % use 16-bit precision to store the complex-valued projections 
            par0.cache_stack_object = par0.save_memory; % store stack_object to disk when no needed 

            
    case 5
        %% lamni chip dataset !! THIS TEST TAKES SEVERAL HOURS and requires full RAM !!
        tested_templates ={'template_tomo_recons_lamino'} ;
%          par0.scanstomo = [984:1850]; %2326];  % smaller angular range 
        par0.scanstomo = [984:3717];   % full angular range 
        par0.tomo_id = []; % Either scan numbers or tomo_id can be given, but not both, if not provided leave tomo_id=[]

        % IO loading 
        par0.fileprefix='';                  % string at the beginning of the filename, related to reconstruction name
        par0.filesuffix = 'test_1';    %% string searched at the end of the data filenames,  No need to add _c or _recons, it will look for it 
        par0.file_extension = 'h5';
        par0.downsample_projections = 0;     % downsample projections by factor of 2^x, set 0 to do nothing and 1,2,.. for different levels of projection binning 
        par0.analysis_path = fullfile(scratch_path, 'data/e17299_lamni_chip_dataset_2018/analysis/');

        par0.clip_amplitude_quantile = 0.95;  %  clip amplitudes in the loaded projections that are exceeding given quantile, if par0.clip_amplitude_quantile == 1, do nothing 
        par0.max_residua_limit = 100;        % limit used to determine which projection have failed 
        par0.surface_calib_file = [];
        par0.omnyposfile = fullfile(scratch_path, 'data/e17299_lamni_chip_dataset_2018/specES1/scan_positions/scan_%05d.dat');		%Filename pattern for Orchestra interferometer position files   
        par0.OMNY_angle_file = fullfile(scratch_path, 'data/e17299_lamni_chip_dataset_2018/specES1/dat-files/tomography_scannumbers.txt');  % Filename with angles

        % Other
        par0.save_memory = true;        % try to limit use of RAM 
            par0.inplace_processing = par0.save_memory; % process object_stack using inplace operations to save memory 
            par0.fp16_precision     = par0.save_memory; % use 16-bit precision to store the complex-valued projections 
            par0.cache_stack_object = par0.save_memory; % store stack_object to disk when no needed 
    otherwise
        error('Missing dataset')
end



for tested_template = tested_templates
    clearvars -except par0 tested_template tested_templates scratch_path GPU_id base_path test_simulated_data test_real_data base_path

    utils.verbose(struct('prefix', 'init'))

    %% test stage 0: load basic configuration parameters 
    utils.verbose(-1,'====================================================')
    utils.verbose(-1,'===== Testing template "%s" ============', tested_template{1})
    utils.verbose(-1,'====================================================')
    
    % set debugging info level and marks 
    debug(1)
    warning('off', 'MATLAB:mpath:nameNonexistentOrNotADirectory')
    warning('off', 'MATLAB:dispatcher:pathWarning')

    run(tested_template{1})

    for item = fieldnames(par0)'
        par.(item{1}) = par0.(item{1}); 
    end
   

    %% test stage 1: load test data and continue with the remplate 
    utils.verbose(-1,'Running template')

    debug(3)
    utils.verbose(struct('prefix', 'template'))
    run(tested_template{1})
    

end

end

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%*-----------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing language this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS matlab package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.

