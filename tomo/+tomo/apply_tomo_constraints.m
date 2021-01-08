% [volume_new, update] = apply_tomo_constraints(volume, mask, lamino_angle , low_freq_protection,  constrain_fun, Niter)
% apply laminography constraints in the real space and try to refill missing cone in laminography by provided prior
%     knowledge 
% Inputs:
%    **volume           - (3D array) represeting the refined volume in realspace 
%    **mask             - (vector, array), mask pushing pixels where mask < 1 towards zero. Can be either 3D or only for example along 3r axis ie size(mask) = [1,1,Nlayers]
%    **lamino_angle     - (scalar), laminography angle from 0 to 90degrees, 90 == classical tomo, it is used to calculate the missing cone  
%    **low_freq_protection - (bool), used to protect in the fourier space the central region, ie low spatial frequncies. Important when multiscale approach is used 
%    **constrain_fun       - anonymous function providing constrains such as positivity or material range limits 
%    **Niter               - number of optimization iterations 
% *returns*
%   ++volume_new             refined object 
%   ++update                 (norm(volume) - norm(update_new)) /   norm(volume)
%
% Example: 
%   see template_tomo_recons_lamino.m for working example 


function [volume_new, update] = apply_tomo_constraints(volume, mask, angles, low_freq_protection,  value_max, value_min, TV_lambda, Niter)
    import utils.Garray 
    
    Npix = size(volume);
    
    fft_mask = tomo.get_tomo_fourier_mask_3d( Npix, angles); 
    fft_mask = Garray(fft_mask);
    
    if low_freq_protection
        % avoid modification of the low spatial frequencies that were
        % already refined 
        fft_mask = fftshift(fft_mask); 
        for i = 1:3
            grid{i} = ceil(Npix(i)/2)+[-ceil(Npix(i)/8):floor(Npix(i)/8)];
        end
        fft_mask(grid{:}) = 0; 
        fft_mask = fftshift(fft_mask); 
    end
    
    volume = Garray(volume); 

    fft_split = 1;
             
    for iter = 1:Niter 
        utils.progressbar(iter,Niter)
        volume_new = volume;

        volume_new = regularization.local_TV3D_chambolle(volume_new, TV_lambda, 10);

        % positivity constraint 
        volume_new = arrayfun(@clip_range,volume_new, value_max, value_min, mask); 
        %volume_new = clip_range(volume_new, value_max, value_min, mask);

        %% go to the Fourier space
        fvolume = (math.fftn_partial(Garray(volume), fft_split));
        fvolume_new = (math.fftn_partial(Garray(volume_new), fft_split));
        
        %% merge updated and original dataset in the fourier space
        %% use overrelaxation of the constraint to get faster convergence 
        relax = 1.5;
        regularize = 0;
        fvolume = arrayfun(@relax_contraint,fvolume, fvolume_new, fft_mask, relax, regularize); 
        clear fvolume_new

        %% back to the real space 
        volume_new = real(math.ifftn_partial(Garray(fvolume), fft_split));
        clear fvolume
        % get difference in update 
        update = gather(norm(volume(:)-volume_new(:)) ./ norm(volume(:)));

        volume = volume_new;  
  
    end
    
    volume = gather(volume); 
    
end

% auxiliary function for fast execution on GPU 
function fvolume = relax_contraint(fvolume, fvolume_new, fft_mask, relax, regularize)
    fvolume = fvolume .* ( 1- relax.*fft_mask) + fvolume_new .* relax.*fft_mask;
    %relax_data = 0.2;
    %fft_mask_data = 1 - fft_mask;
    %fvolume = fvolume .* ( 1- relax_data*fft_mask_data) + fvolume_new .* relax_data.*fft_mask_data;
    fvolume = fvolume .* (1 - regularize.*fft_mask); 
end

function array = clip_range(array, max_val, min_val, mask)
    array = max(min_val, min(max_val, array)) .* mask;
    %array = max(min_val, min(max_val, array));
end

