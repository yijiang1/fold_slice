% Wrapper function for simulating diffraction patterns for electron ptychography 
% Results are saved in formats required by reconstruction scripts
% Assumes the strong phase approximation model only
%
% input:
% ** param structure containting experimental conditions and probe parameters

% returns:
% dp                 simulated diffraction patterns
% probe_true         simulated probe
% p                  p structure corresponding saved probe function

function [dp, probe_true, p] = sim_electron_cbed(param)

    addpath(strcat(pwd,'/utils_electron/'))
    addpath(core.find_base_package)

    % parse inputs
    parse_param = inputParser;
    parse_param.KeepUnmatched = true;
    
    % parameters
    parse_param.addParameter('output_path', '', @ischar)

    parse_param.addParameter('field_of_view', 36, @isnumeric) %scan FOV in angstroms
    parse_param.addParameter('N_dp', 256, @isnumeric) %size of diffraction pattern in pixels 
    parse_param.addParameter('N_dp_orig', 256, @isnumeric) %original size of diffraction pattern in pixels
    parse_param.addParameter('scan_step_size', 2, @isnumeric) %scan step size in angstroms
    parse_param.addParameter('max_position_error', 0, @isnumeric) %max random position error in angstroms
    parse_param.addParameter('dose', 0, @isnumeric) %total electron dose in # of electrons per angstroms^2
    parse_param.addParameter('voltage', 80, @isnumeric)  %beam voltage in keV
    parse_param.addParameter('alpha_max', 20, @isnumeric)  %semi-convergence angle in mrads

    % object parameters
    parse_param.addParameter('object', 0, @isnumeric)  %complex object function
    parse_param.addParameter('dx', 0.1, @isnumeric)  %real-space pixel size in angstroms

    % probe parameters
    parse_param.addParameter('probe_df', 0, @isnumeric) % defocus in  angstroms
    parse_param.addParameter('probe_c3', 0, @isnumeric) % spherical aberrations in angstroms
    parse_param.addParameter('probe_c5', 0, @isnumeric )
    parse_param.addParameter('probe_c7', 0, @isnumeric )
    parse_param.addParameter('probe_f_a2', 0, @isnumeric )
    parse_param.addParameter('probe_theta_a2', 0, @isnumeric )
    parse_param.addParameter('probe_f_a3', 0, @isnumeric )
    parse_param.addParameter('probe_theta_a3', 0, @isnumeric )
    parse_param.addParameter('probe_f_c3', 0, @isnumeric )
    parse_param.addParameter('probe_theta_c3', 0, @isnumeric )
    
    % 
    parse_param.addParameter('gpu_id', 1, @isnumeric )
    parse_param.addParameter('overwrite_data', 0, @islogical )

    % simulation parameters
    %parse_param.addParameter('N_data', 1, @isnumeric) % number of datasets to be simulated

    parse_param.parse(param)
    param_input = parse_param.Results;
    
    %% check inputs
    assert(~isempty(param_input.output_path), 'Output path cannot be empty!')

    assert(param_input.N_dp>0, 'Invalid diffraction pattern size!')
    assert(param_input.field_of_view>0, 'Invalid scan field of view!')
    assert(param_input.scan_step_size>0, 'Invalid scan step size!')
    assert(param_input.voltage>0, 'Invalid voltage!')
    assert(param_input.alpha_max>0, 'Invalid aperature size!')
    assert(param_input.dose>0, 'Invalid dose level!')

    field_of_view = param_input.field_of_view;
    N_dp = param_input.N_dp;
    N_dp_orig = param_input.N_dp_orig;
    scan_step_size = param_input.scan_step_size;
    max_position_error = param_input.max_position_error;
    
    dose = param_input.dose;
    dx = param_input.dx;
    
    %gpuDevice(param_input.gpu_id);
    %% prepare parameters    
    N_scans_h = ceil(field_of_view/scan_step_size); % number of scan positions along horizontal direction
    N_scans_v = ceil(field_of_view/scan_step_size); % number of scan positions along vertical direction

    Nc_avg = dose*scan_step_size^2/N_dp^2; %average electron count per detector pixel. For poisson noise, SNR = sqrt(Nc_avg);

    N_obj = size(param_input.object);
    ind_obj_center = floor(N_obj/2)+1;

    %object = gpuArray(param_input.object);
    object = param_input.object;

    %% generate probe function
    %disp('Generate probe function...')
    par_probe = {};
    par_probe.df = param_input.probe_df; %defocus in angstrom
    par_probe.C3 = param_input.probe_c3; %third-order spherical aberration in angstrom
    par_probe.C5 = param_input.probe_c5; %fifth-order spherical aberration in angstrom
    par_probe.C7 = param_input.probe_c7; %seventh-order spherical aberration in angstrom
    par_probe.f_a2 = param_input.probe_f_a2; %azimuthal orientation in radian
    par_probe.theta_a2 = param_input.probe_theta_a2; %twofold astigmatism in angstrom
    par_probe.f_a3 = param_input.probe_f_a3; %threefold astigmatism in angstrom
    par_probe.theta_a3 = param_input.probe_theta_a3; %azimuthal orientation in radian

    par_probe.voltage = param_input.voltage; %beam voltage in keV
    par_probe.alpha_max = param_input.alpha_max; %semi-convergence angle in mrad
    par_probe.plotting = false;
    dk = 1/dx/N_dp_orig; %fourier-space pixel size in 1/A

    [probe_true, ~] = make_tem_probe(dx, N_dp_orig, par_probe);
    %probe_true = crop_pad(probe_true_init, [N_dp_orig, N_dp_orig]);
	%probe_true = gpuArray(probe_true);

    %% save initial probe and parameters
    probe = crop_pad(probe_true, [N_dp, N_dp]);
    if ~exist(param_input.output_path, 'dir'); mkdir(param_input.output_path); end
    
    p = {};
    p.binning = false;
    p.detector.binning = false;
    p.dk = dk*N_dp_orig/N_dp;
    p.N_scans_h = N_scans_h;
    p.N_scans_v = N_scans_v;
    save(fullfile(param_input.output_path, 'init_probe'), 'probe', 'p', 'par_probe')
    save(fullfile(param_input.output_path, 'true_probe'), 'probe_true', 'par_probe')

    %% Generate scan positions
    pos_h = (1+(0:N_scans_h-1)*param_input.scan_step_size);
    pos_v = (1+(0:N_scans_v-1)*param_input.scan_step_size);
    % centre this
    pos_h  = pos_h - (mean(pos_h));
    pos_v  = pos_v - (mean(pos_v));
    [Y,X] = meshgrid(pos_h, pos_v);
    Y = Y';
    X = X';
    pos_true_h = X(:); % true posoition
    pos_true_v = Y(:);

    %pos_recon_init_h = pos_true_h;
    %pos_recon_init_v = pos_true_v;

    %add random position errors - to simulate scan noise
    pos_true_h = pos_true_h + max_position_error*(rand(size(pos_true_h))*2-1);
    pos_true_v = pos_true_v + max_position_error*(rand(size(pos_true_v))*2-1);

    %calculate indicies for all scans
    N_scan = length(pos_true_h);
    %position = pi(integer) + pf(fraction)
    pv_i = round(pos_true_v/dx);
    pv_f = pos_true_v - pv_i*dx;
    ph_i = round(pos_true_h/dx);
    ph_f = pos_true_h - ph_i*dx;

    ind_h_lb = ph_i - floor(N_dp_orig/2) + ind_obj_center(2);
    ind_h_ub = ph_i + ceil(N_dp_orig/2) -1 + ind_obj_center(2);
    ind_v_lb = pv_i - floor(N_dp_orig/2) + ind_obj_center(1);
    ind_v_ub = pv_i + ceil(N_dp_orig/2) -1 + ind_obj_center(1);
    
    %% generate diffraction patterns
    %dp = gpuArray(zeros(N_dp, N_dp, N_scan));
    dp = zeros(N_dp, N_dp, N_scan);
    
    if isfile(fullfile(param_input.output_path, 'data_roi0_dp.hdf5')) && param_input.overwrite_data
        delete(fullfile(param_input.output_path, 'data_roi0_dp.hdf5'));
        delete(fullfile(param_input.output_path, 'data_roi0_para.hdf5'));
    elseif isfile(fullfile(param_input.output_path, 'data_roi0_dp.hdf5'))
        return
    end
    
    %snr = ones(N_scan, 1)*inf; %signal-to-noise ratio of each diffraction pattern
    for i=1:N_scan
        probe_s = shift(probe_true, dx, dx, ph_f(i), pv_f(i));
        object_roi = object(ind_v_lb(i):ind_v_ub(i),ind_h_lb(i):ind_h_ub(i));
        psi =  object_roi .* probe_s;

        %FFT to get diffraction pattern
        dp_true = abs(fftshift(fft2(ifftshift(psi)))).^2;

        %dp_true = imresize(dp_true, N_dp/N_dp_orig, 'box');
        dp_true = dp_true(1:N_dp_orig/N_dp:end, 1:N_dp_orig/N_dp:end);
        %dp_true(dp_true<0) = 0;
        dp_temp = dp_true;
        
        %Add poisson noise
        if Nc_avg<inf
            dp_temp = dp_true/sum(dp_true(:))*(N_dp^2*Nc_avg);
            
            dp_temp = poissrnd(dp_temp);
            
            dp_temp = dp_temp*sum(dp_true(:))/(N_dp^2*Nc_avg);
            %snr(i) = mean((dp_true(:)))/std(dp_true(:) - dp_temp(:));
        end
        
        dp(:,:,i) = dp_temp;
    end
    
    %dp = gather(dp);

    % save diffraction patterns 
    %save_dir = fullfile(base_path, data_path, strcat('data',num2str(j)));
    if ~exist(param_input.output_path, 'dir'); mkdir(param_input.output_path); end
    
    h5create(fullfile(param_input.output_path, 'data_roi0_dp.hdf5'), '/dp', size(dp), 'ChunkSize',[size(dp,1) size(dp,1), 1], 'Deflate',4)
    h5write(fullfile(param_input.output_path, 'data_roi0_dp.hdf5'), '/dp', dp)

    h5create(fullfile(param_input.output_path, 'data_roi0_para.hdf5'), '/ppX', size(pos_true_h))
    h5write(fullfile(param_input.output_path, 'data_roi0_para.hdf5'), '/ppX', pos_true_h)
    h5create(fullfile(param_input.output_path, 'data_roi0_para.hdf5'), '/ppY', size(pos_true_v))
    h5write(fullfile(param_input.output_path, 'data_roi0_para.hdf5'), '/ppY', pos_true_v)
    
end
