% Wrapper function for reconstructing experimental ptychography data
% Can be used with BO-GP to calibrate experimental parameters
% 
% Output:
%  ++ data_error - averaged data error of last iteration

function [data_error] = ptycho_recon_pear(params, varargin)
    
    parser = inputParser;

    parser.addParameter('defocus',  0 , @isnumeric)  %STEM probe's defocus in angstroms
    parser.addParameter('thickness',  0 , @isnumeric)  %total sample thickness in angstroms for electron or meters for x-ray
    parser.addParameter('rot_ang',  0 , @isnumeric)  %in-plane rotation between diffraction pattern and scanning plane in degrees
    parser.addParameter('alpha_max',  0 , @isnumeric)  %aperture size in mrad
    parser.addParameter('beam_source', '', @ischar) %for e-ptycho.
    parser.addParameter('zp_propagation_distance',  0 , @isnumeric)  %x-ray probe's propagation distance in meters
    
    parser.addParameter('gpu_id',  1 , @isnumeric)  %GPU ID
    parser.addParameter('use_external_probe', 0, @islogical)
    parser.addParameter('initial_probe_file', '', @ischar)
    parser.addParameter('use_external_object', 0, @islogical)
    parser.addParameter('initial_object_file', '', @ischar)
    parser.addParameter('use_external_positions', 0, @islogical)
    parser.addParameter('initial_position_file', '', @ischar)
    parser.addParameter('use_grid_scan', false, @islogical)

    parser.addParameter('multislice_ptycho', 0, @islogical)
    parser.addParameter('position_corr', 0, @islogical)
    parser.addParameter('variable_probe_corr', 0, @islogical)
    parser.addParameter('variable_intensity_corr', 0, @islogical)
    parser.addParameter('multimodal_update', 0, @islogical)
    parser.addParameter('use_momentum', 0, @islogical)

    parser.parse(varargin{:})
    r = parser.Results;

    % load all to the param structure
    par = params;
    for name = fieldnames(r)'
        if ~isfield(par, name{1}) || ~ismember(name, parser.UsingDefaults) % prefer values in param structure 
            par.(name{1}) = r.(name{1});
        end
    end
    
    par_recon = par;
    par_recon.gpu_id = par.gpu_id; 
    
    if strcmp(par.beam_source, 'electron') && ~isfield(par, 'dk')
        par_recon.d_alpha = par.alpha_max / par.rbf; %use rbf to calculate pixel size of diffraction patterns
    end
    
    affine_mat  = compose_affine_matrix(1, 0, par.rot_ang, 0);
    par_recon.affine_matrix_11 = affine_mat(1,1);
    par_recon.affine_matrix_12 = affine_mat(1,2);
    par_recon.affine_matrix_21 = affine_mat(2,1);
    par_recon.affine_matrix_22 = affine_mat(2,2);

    % initial probe
    if par_recon.use_external_probe
        par_recon.use_model_probe = false;
        par_recon.initial_probe_file = par.initial_probe_file;
        par_recon.normalize_init_probe = false;
    else
        par_recon.use_model_probe = true;
        par_recon.initial_probe_file = '';
        par_recon.normalize_init_probe = true;
        %xray
        par_recon.model_probe_prop_dist = par.zp_propagation_distance;
    end
    
    %initial object
    if par_recon.use_external_object
        par_recon.use_model_object = false;
        par_recon.initial_object_file = par.initial_object_file;
        par_recon.multiple_layers_obj = true;
    else
        par_recon.use_model_object = true;
        par_recon.initial_object_file = '';
    end
    
    %initial positions
    if par_recon.use_external_positions
        par_recon.src_positions = 'matlab_pos';
        par_recon.scan_type = 'custom_GPU';
        par_recon.custom_positions_source = par.initial_position_file;
    else
        if par_recon.use_grid_scan
            par_recon.src_positions = 'matlab_pos';
            par_recon.scan_type = 'raster';
        else
            par_recon.src_positions = 'hdf5_pos'; % Load scan positions from an hdf5 file
            par_recon.scan_type = 'default'; % Type of scan positions
        end
    end
    
    if par.position_corr; par_recon.probe_position_search = 0; end
    if par.variable_probe_corr; par_recon.variable_probe_modes = 1; end
    if par.variable_intensity_corr; par_recon.variable_intensity = true; end
    if par.multimodal_update; par_recon.apply_multimodal_update = true; end
    if par.use_momentum
        par_recon.method = 'MLc';
        par_recon.momentum = 0.5;
    end
    
    par_recon.probe_alpha_max = par.alpha_max;
    par_recon.probe_df = par.defocus;

    par_recon.output_dir_suffix = '_pear';
    
    if par.multislice_ptycho
        par_recon.eng_name = 'GPU_MS';
        par_recon.delta_z = par.thickness / par_recon.Nlayers;
    end
    
    [~, ~, data_error] = ptycho_recon(par_recon);
    
end
