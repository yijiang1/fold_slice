function [sinogram] = get_sinogram_from_complex_object(object, par)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
    
    
    if isfield(par,'par.complex_proj_close') && ~isempty(par.complex_proj_close)
        object_mask = tomo.estimate_reliability_region_grad(object, par.complex_proj_close, par.complex_proj_erode, par.proj_mask_par);
    else
        object_mask = 1;
    end
    
    sinogram = -tomo.unwrap2D_fft2_split(object, [], 1, object_mask, par.GPU_list);
    
    % apply a mask to phase
    if par.get_sino_weights_from_object
        sinogram_msak = tomo.estimate_reliability_region_grad(object, par.real_proj_close, par.real_proj_erode, par.proj_mask_par);
    else
        sinogram_msak = tomo.estimate_reliability_region_grad(sinogram, par.real_proj_close, par.real_proj_erode, par.proj_mask_par);
    end

    % apply mask to sinogram
    %sino_weights = sino_weights(selected_ROI{1},selected_ROI{2},:);
    sinogram = sinogram .* sinogram_msak;

end


