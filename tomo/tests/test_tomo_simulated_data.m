%% TEMPLATE FOR AUTOMATIC TOMOGRAPHY CODE TESTS 
% test uses synthetic data to run all templates under well controlled
% conditions, use setting if you want to add noise or other difficulties to
% be tested during alignment 
cd(fullfile( fileparts(mfilename('fullpath')), '..'))
addpath('tests')
addpath('utils')
addpath('./')
addpath(find_base_package)
 
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Edit this section %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

tested_templates = 1:3;  % IDs of the tested templates 

templates_list ={   'template_tomo_recons', ...
                    'template_tomo_recons_lamino', ...
                    'template_tomo_nonrigid', ...
                    'template_tomo_recons_deprecated'...
                  };

      
verbose_level = -1;    % -1 = keep very quiet the reconstructions 

% artificial data settings
Npix_vol = [200,200,200] ;
undersampling = 1;              % level of undersampling compared to Crowther criterion 
noise_level = 0;                % relative noise level with respect to the maximal phase value in the projections 
N_subtomos = 4; 
add_residual_layer = false;     % add a thin metal-like layer to test behaviour with residua
Nangles = ceil(pi/2*Npix_vol(1)) / undersampling; 
Nangles = ceil(Nangles/N_subtomos)*N_subtomos; % make splitable for 4 subtomos
asize = ceil(Npix_vol(1:2) / 10);  % size of the simulated probe 

if ~exist('GPU_id', 'var'); GPU_id = [1]; end
if ~exist('base_path', 'var'); base_path = '../'; end
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%



assert(all(GPU_id <= gpuDeviceCount), 'Select valid GPU id in GPU_id')

utils.verbose(verbose_level)
setenv('TMP',[base_path,'/tmp']) % TEMP for matlab scripts 



%% MAKE A PHANTOM
utils.verbose(-1,'Creating phantom')
rng default
volData = randn(Npix_vol - asize(1), 'single'); 
volData = (0.5+0.5*(utils.imgaussfilt3_fft(volData, 3) > 0)) .* ...
          (utils.imgaussfilt3_fft(volData, 10) > 0);
volData = utils.imgaussfilt3_conv(volData, 0.6); % prevent too sharp edges 
      
if add_residual_layer
    % add kind of metalic layer to test robustness of the codes 
    [X,Y,Z] = meshgrid(-Npix_vol(1)/2:Npix_vol(1)/2-1,-Npix_vol(2)/2:Npix_vol(2)/2-1,-Npix_vol(3)/2:Npix_vol(3)/2-1); 
    layer = abs(0.1*X + 0.5*Y+0.2*Z + 50*utils.imgaussfilt3_fft(randn(Npix_vol, 'single'),10) )<0.5; 
    volData = volData + 10*utils.crop_pad_3D(layer, size(volData)); 
end
volData = utils.apply_3D_apodization(volData, 0, Npix_vol(1)/5, 1); 
volData = utils.crop_pad_3D(volData, Npix_vol);

for tested_template = templates_list(tested_templates)
    utils.verbose(struct('prefix', 'init'))

    %% test stage 0: load basic configuration parameters 
    utils.verbose(-1,'Loading template "%s"', tested_template{1})
    if strcmpi(tested_template{1}, 'template_tomo_nonrigid')
        %nonrigid tomo has data generation already included in the template 
        debug(1)
        run(tested_template{1})
    else
        par = struct();
        par.lamino_angle = 90;  % default setting 

        % set debugging info level and marks 
        debug(1)
        warning('off', 'MATLAB:mpath:nameNonexistentOrNotADirectory')
        warning('off', 'MATLAB:dispatcher:pathWarning')

        run(tested_template{1})

        utils.verbose(-1,'Creating data')

        % create "data"
        if ~par.is_laminography
            theta = pi+linspace(0, 180*(1-1/Nangles), Nangles); 
        else
            Nangles = Nangles * 2;
            theta = pi+linspace(0, 360*(1-1/Nangles), Nangles); 
        end

        % create "N_subtomos" subtomos 
        theta = reshape(reshape(theta, N_subtomos,[])',1,[]); 
        par.subtomos = reshape(ones(Nangles/N_subtomos,N_subtomos).*[1:N_subtomos],1,[]); 
        
        % Add noise and offsets to the measured positions to make the
        % alignment more difficult 
        position_errors = (0.1*randn(Nangles,2) + 0.1*sind(3*theta')) * Npix_vol(3); 
 
        asize = asize + ceil(max(asize, max(abs(position_errors)))/32)*32;  % avoid the sample going out of FOV
        par.asize = asize; 

        if par.lamino_angle == 90
            Npix_proj = [Npix_vol(3), ceil(sqrt(2)*Npix_vol(1))]+asize; 
            max_sample_height = inf; 
        else
            Npix_proj = ceil(sqrt(2)*Npix_vol([3,1]) .* [cosd(par.lamino_angle), 1] )+asize; 
            max_sample_height = 5e-6 ; 
        end



 
        % provide other parameters required by the template 
        par.factor = 1; 
        par.factor_edensity = 1;
        par.pixel_size = 50e-9;
        par.scans_string = '';
        par.output_folder = '???';  % some nonexistent folder, in default the template should not write any data during debug mode 
        obj_interf_pos_x = 0; 
        obj_interf_pos_y = 0; 
        delta_stack_prealign = [] ; 
        par.tilt_angle = 0;
        par.skewness_angle = 0; 
        par.lambda = 0.2e-9; 
        par.scanstomo = 1:Nangles;
        par.air_gap = [5,5];
        par.GPU_list = GPU_id; 
        par.nresidua_per_frame = 0;
        
        % generate geometry 
        Npix_vol(3) = min(Npix_vol(3),ceil(max_sample_height / par.pixel_size)-1); 
        
        CoR = Npix_proj/2 + position_errors; 
        [cfg, vectors]  = astra.ASTRA_initialize(Npix_vol, Npix_proj, theta, par.lamino_angle, 'rotation_center', CoR); 
        % find optimal split, for small volumes below 600^3 no split is needed 
        split  = astra.ASTRA_find_optimal_split(cfg); 
        par.illum_sum = ones(Npix_proj);

        
        % generate complex projections 
        stack_object  = tomo.Ax_sup_partial(utils.crop_pad_3D(volData, Npix_vol), cfg, vectors, split); 
        stack_object = stack_object / math.sp_quantile(stack_object, 0.99, 10); 
        if noise_level>0; stack_object = stack_object + noise_level * randn(size(stack_object)); end
        stack_object = exp(-0.1*stack_object - 4i*stack_object);

        object = stack_object(:,:,1);

        
        if strcmpi(tested_template, 'template_tomo_recons_lamino')
            % create errors in global geometry -> test automatic refinement
            par.tilt_angle = -0.5;
            par.skewness_angle = 0.5; 
        else
            par.tilt_angle = 0;
            par.skewness_angle = 0; 
        end


        %% test stage 1: load test data and continue with the remplate 
        utils.verbose(-1,'Loading template "%s"', tested_template{1})

        debug(2)
        utils.verbose(struct('prefix', 'template'))
        run(tested_template{1})
        
        
        %   check results of the geometry refinement provided by tests 
        if strcmpi(tested_template, 'template_tomo_recons_lamino')
            assert(abs(par.tilt_angle) < 0.1 && abs(par.skewness_angle) < 0.1, 'Geometry refinement in laminography did not work well')
        end
    end
    clearvars -except par tested_templates Nangles volData Npix_vol asize N_subtomos GPU_id noise_level test_astra test_mex_function test_simulated_data test_real_data base_path 


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

