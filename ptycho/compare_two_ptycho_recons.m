% Wrapper function for comparing two ptychographic reconstructions
% Can be called in python and used for automatic parameter tuning
%
% Input:
% ** param    - structure containting parameters for reconstruction files and
% evaluation method
%
% Outputs:
% ++ score    - similarity score between the input reconstructions

function [score] = compare_two_ptycho_recons(param)
    
    addpath ../
    addpath('utils')

    params = struct;
    % parse inputs
    parse_param = inputParser;
    
    parse_param.addParameter('base_path',  '' , @ischar)
    parse_param.addParameter('scan_string_format',  '' , @ischar)
    parse_param.addParameter('roi_label',  '' , @ischar)
    parse_param.addParameter('recon_dir',  '' , @ischar)
    parse_param.addParameter('verbose_level',  0 , @isnumeric)
    parse_param.addParameter('scan1', 1, @isnumeric)
    parse_param.addParameter('scan2', 2, @isnumeric)

    parse_param.addParameter('thickring', 3, @isnumeric)
    parse_param.addParameter('taper', 20, @isnumeric)
    parse_param.addParameter('SNRt', 0.5, @isnumeric)
    
    parse_param.addParameter('crop_roi_y_lb', -1, @isnumeric)
    parse_param.addParameter('crop_roi_y_ub', -1, @isnumeric)
    parse_param.addParameter('crop_roi_x_lb', -1, @isnumeric)
    parse_param.addParameter('crop_roi_x_ub', -1, @isnumeric)
    
    parse_param.addParameter('offset_x', 0, @isnumeric)
    parse_param.addParameter('offset_y', 0, @isnumeric)
    
    parse_param.addParameter('crop_x', 0, @isnumeric)
    parse_param.addParameter('crop_y', 0, @isnumeric)
    
    parse_param.addParameter('electron', 0, @islogical)

    parse_param.addParameter('metric', 'ssim', @ischar)
    parse_param.addParameter('file1', '', @ischar)
    parse_param.addParameter('file2', '', @ischar)

    parse_param.addParameter('lambda_phase_range', 0, @isnumeric)

    parse_param.parse(param)
    param_input = parse_param.Results;

    %% check inputs
    %assert(~isempty(param_input.base_path), 'Base path cannot be empty!')
    %assert(~isempty(param_input.recon_dir), 'Recon path cannot be empty!')

    %% For .mat reconstruction output
    base_path = param_input.base_path;
    recon_dir = param_input.recon_dir;
    roi_label = param_input.roi_label;
    scan_string_format = param_input.scan_string_format;
    scanNo_s = [param_input.scan1, param_input.scan2];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Alignment parameters %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    params.verbose_level = param_input.verbose_level;   % adjust output level
    params.plotting = params.verbose_level;        % (3) show everything, (2) show aligned images + FSC, (1) show FSC, (0) none
    params.remove_ramp = 1;     % Try to remove ramp from whole image before initial alignment
    params.image_prop = 'phasor'; % = 'complex' or = 'phasor' (phase with unit amplitude) or = 'phase'  (Note: phase should not be used if there is phase wrapping)
    params.crop = ''; 
                % '' for using the default half size of the probe
                % 'manual' for using GUI to select region. This will display the range, e.g. {600:800, 600:800}     
                % {600:800, 600:800} for custom vertical and horizontal cropping, respectively
    params.flipped_images = 0; % If images are taken with a horizontal flip, e.g. 0 & 180 for tomography
    params.GUIguess = 0;       % To click for an initial alignment guess, ignores the values below
    params.guessx = [];       % Some initial guess for x alignment
    params.guessy = [];
    params.electron = param_input.electron;

    %%%%%%%%%%%%%%%%%%%%%%
    %%% FSC parameters %%%
    %%%%%%%%%%%%%%%%%%%%%%
    params.taper = param_input.taper;             % Pixels of image tapering (smoothing at edges) - Increase until the FSC does not change anymore
    params.SNRt = param_input.SNRt;      % SNRt = 0.2071 for 1/2 bit threshold for resolution of the average of the 2 images
                                         % SNRt = 0.5    for 1   bit threshold for resolution of each individual image
    params.thickring = param_input.thickring;  % Thickness of Fourier domain ring for FSC in pixels
    params.freq_thr = 0.05;  % (default 0.05) To ignore the crossings before freq_thr for determining resolution 

    %%%%%%%%%%%%
    %%% misc %%%
    %%%%%%%%%%%%
    params.prop_obj = false;     % propagation distance at the sample plane; leave empty to use the value from the reconstruction p structure; set to "false" for no propagation
    params.apod = [];            % if true, applies an apodization before propagating by params.prop_obj, the apodization border is around the valid reconstruction region; leave empty to use the value from the reconstruction p structure
    params.lambda = [];           % wavelength; needed for propagating the object; leave empty to use the value from the reconstruction p structure
    params.pixel_size = [];       % pixel size at the object plane; leave empty to use the value from the reconstruction p structure

    %%%%%%%%%%%%%%%%%%%%
    %%% FP parameter %%%
    %%%%%%%%%%%%%%%%%%%%

    %%% the following parameters are ignored, unless p.fourier_ptycho==true %%%
    params.filter_FFT = true;           % apply a circular mask to the reconstructed spectrum (needs p.plot_maskdim)
    params.crop_factor = 0.9;           % crop final image by the given factor
    params.crop_asize = [800 800];      % crop object before applying the FFT
    params.z_lens = 49.456e-3;          % FZP focal distance

    %% For .mat file
    file = {};
    for i=1:2
        if i==1 && ~isempty(param_input.file1)
            recon_file = param_input.file1;
        elseif i==2 && ~isempty(param_input.file2)
            recon_file = param_input.file2;
        else
            sub_dir = strcat(sprintf(scan_string_format, scanNo_s(i)), '/roi' ,roi_label, '/');
            recon_file = fullfile(base_path, sub_dir, recon_dir, 'Niter*');
        end
        recon_file = dir(recon_file);

        if ~isempty(recon_file)
            recon_folder = recon_file(end).folder;
            recon_file = recon_file(end).name;
            recon_file = fullfile(recon_folder, recon_file);
            load(recon_file)
            file{i} = object;
            params.pixel_size = p.dx_spec;
            if param_input.crop_roi_y_lb>0 || param_input.crop_roi_x_lb>0
                crop_roi{1} = param_input.crop_roi_y_lb:param_input.crop_roi_y_ub;
                crop_roi{2} = param_input.crop_roi_x_lb:param_input.crop_roi_x_ub;
            else
                if param_input.crop_y~=0 || param_input.crop_x~=0
                    cen = floor(size(object)/2)+1;
                    crop_roi{1} = cen(1)-param_input.crop_y:cen(1)+param_input.crop_y-1;
                    crop_roi{2} = cen(2)-param_input.crop_x:cen(2)+param_input.crop_x-1;
                else
                    offset_y = param_input.offset_y;
                    offset_x = param_input.offset_x;

                    crop_roi{1} = p.object_ROI{1}(1)+offset_y:p.object_ROI{1}(end)-offset_y;
                    crop_roi{2} = p.object_ROI{2}(1)+offset_x:p.object_ROI{2}(end)-offset_x;
                end
            end
            params.crop{i} = crop_roi;       
        else
            disp(fullfile(base_path,sub_dir,recon_dir))
            error('No recon file found')
        end
    end
    [resolution, stat, subim1, subim2] = aligned_FSC(file{1}, file{2}, params);
    
    subim1_ph = phase_unwrap(angle(subim1));
    subim2_ph = phase_unwrap(angle(subim2));
	
    switch param_input.metric
        case 'frc_1bit'
            score = resolution(1);
        case 'frc_AUC'
            score = 1 - stat.fsc_mean;
        case 'ssim'
            score = ssim(subim1_ph, subim2_ph);
        case 'psnr'
            score = psnr(subim1_ph, subim2_ph);
    end
    
    if param_input.lambda_phase_range~=0
        subim1_max = max(angle(subim1(:)));
        subim1_min = min(angle(subim1(:)));
        subim1_phase_range = subim1_max-subim1_min;
        subim2_max = max(angle(subim2(:)));
        subim2_min = min(angle(subim2(:)));
        subim2_phase_range = subim2_max-subim2_min;
        subims_phase_range = (subim1_phase_range + subim2_phase_range)/2;
        score = score + subims_phase_range * param_input.lambda_phase_range;
    end
end
