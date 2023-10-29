% Wrapper function for reconstructing experimental electron ptychography data
% Can be used with BO-GP to tune experimental parameters
% Written by Yi Jiang
%
% Inputs:
% ** par       - structure containting data and reonstruction parameters
% ** defocus   - probe defocus (in angstroms)
% ** alpha     - aperture size (in mrad)
% ** rot_ang   - in-plane rotation between diffraction patterns and scan cor
% ** thickness - total sample thickness (in angstroms). Effective in multislice recon only
% ** GPU_list  - list of GPU ids. Assumes multiple workers if length(GPU_list) > 1
%
% Output:
%  ++ data_error - averaged data error of last iteration

function [data_error] = ptycho_recon_exp_data_electron(par, defocus, alpha, rot_ang, thickness, GPU_list)

if length(GPU_list)>1
    t = getCurrentTask;
    t = t.ID;
    gpu_id = GPU_list(t);
else
    gpu_id = GPU_list;
end

par_recon = par;
par_recon.gpu_id = gpu_id; 

if ~isfield(par, 'dk') || par_recon.dk<0
    par_recon.d_alpha = alpha/par.rbf; %use rbf to calculate pixel size of diffraction patterns
end

affine_mat  = compose_affine_matrix(1, 0, rot_ang, 0);
par_recon.affine_matrix_11 = affine_mat(1,1);
par_recon.affine_matrix_12 = affine_mat(1,2);
par_recon.affine_matrix_21 = affine_mat(2,1);
par_recon.affine_matrix_22 = affine_mat(2,2);

if par_recon.use_model_probe
    par_recon.probe_alpha_max = alpha;
    par_recon.probe_df = defocus;
end

par_recon.output_dir_suffix = sprintf('_alpha%0.2f_df%0.2f_rot_ang%0.1f', alpha, defocus, rot_ang);

if strcmp(par_recon.eng_name, 'GPU_MS') && thickness > 0
    par_recon.delta_z = thickness / par_recon.Nlayers;
    par_recon.output_dir_suffix = strcat(par_recon.output_dir_suffix, sprintf('_thickness%0.1f', thickness));
end

[~, ~, data_error] = ptycho_recon(par_recon);

end


