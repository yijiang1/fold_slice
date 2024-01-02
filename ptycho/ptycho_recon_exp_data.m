% Wrapper function for reconstructing experimental ptychography data
% Can be used with BO-GP to calibrate experimental parameters
% 
% Output:
%  ++ data_error - averaged data error of last iteration

function [data_error] = ptycho_recon_exp_data(params, varargin)
    
    parser = inputParser;

    parser.addParameter('defocus',  0 , @isnumeric)  %STEM probe's defocus in angstroms
    parser.addParameter('thickness',  0 , @isnumeric)  %total sample thickness in angstroms for electron or meters for x-ray
    parser.addParameter('rot_ang',  0 , @isnumeric)  %in-plane rotation between diffraction pattern and scanning plane in degrees
    parser.addParameter('alpha_max',  0 , @isnumeric)  %aperture size in mrad
    parser.addParameter('beam_source', '', @ischar) %for e-ptycho.
    parser.addParameter('probe_prop_dist',  0 , @isnumeric)  %x-ray probe's propagation distance in meters

    parser.parse(varargin{:})
    r = parser.Results;

    % load all to the param structure
    par = params;
    for name = fieldnames(r)'
        if ~isfield(par, name{1}) || ~ismember(name, parser.UsingDefaults) % prefer values in param structure 
            par.(name{1}) = r.(name{1});
        end
    end

    if length(par.GPU_list)>1 %assume parallel processing
        t = getCurrentTask;
        t = t.ID;
        gpu_id = par.GPU_list(t);
    else
        gpu_id = par.GPU_list;
    end
    
    par_recon = par;
    par_recon.gpu_id = gpu_id; 
    
    if strcmp(par.beam_source, 'electron') && ~isfield(par, 'dk')
        par_recon.d_alpha = par.alpha_max / par.rbf; %use rbf to calculate pixel size of diffraction patterns
    end
    
    affine_mat  = compose_affine_matrix(1, 0, par.rot_ang, 0);
    par_recon.affine_matrix_11 = affine_mat(1,1);
    par_recon.affine_matrix_12 = affine_mat(1,2);
    par_recon.affine_matrix_21 = affine_mat(2,1);
    par_recon.affine_matrix_22 = affine_mat(2,2);
    
    par_recon.probe_alpha_max = par.alpha_max;
    par_recon.probe_df = par.defocus;
    par_recon.model_probe_prop_dist = par.probe_prop_dist;

    par_recon.output_dir_suffix = generate_output_dir_suffix(par.output_dir_suffix_base, varargin, strcmp(par.beam_source, 'electron'));
    
    if strcmp(par_recon.eng_name, 'GPU_MS') && par.thickness > 0
        par_recon.delta_z = par.thickness / par_recon.Nlayers;
    end
    
    [~, ~, data_error] = ptycho_recon(par_recon);

end


function [output_dir_suffix] = generate_output_dir_suffix(output_dir_suffix, param_var, electron)

    if electron %electron ptycho
        for i=1:2:length(param_var)
            switch param_var{i}
                case 'defocus'
                    par_format = '_df%0.2fA';
                case 'thickness'
                    par_format = '_thickness%0.1fA';
                case 'rot_ang'
                    par_format = '_rot_ang%0.1f';
                case 'alpha_max'
                    par_format = '_alpha%0.2mrad';
            end
            output_dir_suffix = sprintf([output_dir_suffix, par_format], param_var{i+1});
        end

    else %X-ray ptycho
        for i=1:2:length(param_var)
            switch param_var{i}
                case 'probe_prop_dist'
                    par_format = '_prop_dist%0.2fm';
                case 'thickness'
                    par_format = '_thickness%0.1fm';
                case 'rot_ang'
                    par_format = '_rot_ang%0.1f';
            end
            output_dir_suffix = sprintf([output_dir_suffix, par_format], param_var{i+1});
        end
        
    end

end
