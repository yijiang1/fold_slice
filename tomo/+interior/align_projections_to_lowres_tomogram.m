% ALIGN_PROJECTIONS_TO_LOWRES_TOMOGRAM align projection from e.g. local tomography to low resolution
%
%  [phase_diff_interior, phase_diff_lres_model, weight, tomogram_lowres] = align_projections_to_lowres_tomogram(stack_object, lres_tomo_path, par)
%
% estimate alignement for the interior tomogram. Resulting alignment is close to optimal
% but it should be further refined using an self-consistent method 
%
% Inputs:
%      **stack_object   - complex valued projection from ptychography 
%      **lres_tomo_path - (string) path to nearfield low resolution
%                   tomogram saved as a .mat file (saved by function tomo.save_tomogram)
%                   The low resolution tomogram assumed to be saved as delta 
%                   (real part of refractive index) 
%                   and it has to include "par" structure with pixel scale
%                   and a conversion factor from delta to phase
%      **par -      tomography parameter structure
%  
% *returns*:
%       ++stack_object - aligned complex valued projections of the high resolution sinogram 
%       ++phase_diff_lres_model - phase difference fot the low resolution preview 
%       ++weight - weights between phase_diff_interior and phase_diff_lres_model, weights are provided as uint8 to save memory 
%       ++shifts - shifts that need to be applied on the phase_diff_interior to be aligned with the phase_diff_lres_model
%       ++resolution_ratio - low resolution pixel size divided by the interior pixel size 
%       ++tomogram_lowres - low resolution tomogram


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



function [stack_object, tomo_lres, win_small, shift, resolution_ratio, Nw_lres] = ...
                    align_projections_to_lowres_tomogram(stack_object, lres_tomo_path, theta, par)

    import utils.* 
    import math.* 
    utils.verbose(struct('prefix', 'prealign'))


    gpu  = gpuDevice();
    if ~isempty(par.GPU_list) && gpu.Index ~= par.GPU_list(1)
        % switch and !! reset !! GPU 
        gpu  = gpuDevice(par.GPU_list(1));
    end
    
    
    %% load low resolution tomogram 
    d  = load(lres_tomo_path);
    tomo_lres = -d.tomogram_delta ./ d.par.factor ;  % revert delta back to phase 
    T = -graythresh(-tomo_lres(tomo_lres<0));
    % calculate the maximal diamter of the local tomogram;
    D = max(sum(radon(max(tomo_lres < T,[],3), 0:180)>0));
    % calculate resolution ratio between low/ high res tomogram 
    resolution_ratio = d.par.pixel_size / par.pixel_size;
    clear d
    
    

    [Nx,Ny,Nangles] = size(stack_object);
    Npix_lres = size(tomo_lres);
    
    % size of the low resolution projections to be generated 
    Nw_lres = ceil([Npix_lres(3), 1.1*D]); % add 10% extra to the maximal diameter 
    Nw_full = ceil(Nw_lres*resolution_ratio); % get the corresponding size of the low res tomogram would be measured in full resolution 

    
    
    %% %%%%%%%%  Get computed sinogram 
    
    % get rough and rather emptirical estimation of the reliability of the interior tomograms 
    win = Garray(single(par.illum_sum)); 
    win = utils.imgaussfilt2_fft(sqrt(win), par.asize(1)/20);
    win = 2-2./(1+ win.^2 / max(win(:).^2));  % limit the weights to 0-1 range 
    win = win .* tukeywin(Nx) .*  tukeywin(Ny)';

   
    gtomo_lres = Garray(tomo_lres);
    
    %% center properly the reconstruction 
    for ii = 1:5
        [x,y,mass] = center(sqrt(max(0,-gtomo_lres))+eps); % abs seems to be more stable than max(0,x) even for missing wedge or laminography 
        % more robust estimation of center 
        rec_center(1) = gather(mean(x.*mass)./mean(mass)); 
        rec_center(2) = gather(mean(y.*mass)./mean(mass));

        % avoid drifts of the reconstructed volume 
        gtomo_lres = tomo.block_fun(@imshift_fft, gtomo_lres, -rec_center(1), -rec_center(2)); 
    end
    tomo_lres = gather(gtomo_lres);
    
    
    verbose(-1,'Aligning projections to a low resolution tomogram')
    
    % get block size to work roughly with 0.2GB arrays 
    Nblocks = ceil((prod(Nw_full)*Nangles*4*2*12)/gpu.AvailableMemory);

    % prepare the local tomo object downsampled to resolution of the low resolution tomogram 
    win_small = interpolate_linear(win,ceil([Nx,Ny]/resolution_ratio));
    win_small = gather(uint8(win_small*255)); 
    
    [stack_object, shift] = ...
        tomo.block_fun(@find_alignment,stack_object, gtomo_lres,theta', par, resolution_ratio, win_small, Nw_lres, struct('Nblocks', Nblocks));
    
    
    verbose(-1,'Pre-alignment done')

    shift = gather(shift);
    
    %% REPORT ALIGNMENT RESULTS 
    figure()
    range = [round(Ny-par.asize(2))/2, round(Nx-par.asize(1))/2]; 
    range = repmat(range, Nangles,1); 
    subplot(1,2,1)
    errorbar(theta, shift(:,1)*par.pixel_size*1e6, range(:,1)*par.pixel_size*1e6 , '.')
    title('Horizontal shift and FOV')
    axis tight
    ylabel('Estimated shift [um]')
    xlabel('Angle [deg]')
    grid on 

    subplot(1,2,2)   
    errorbar(theta, shift(:,2)*par.pixel_size*1e6, range(:,2)*par.pixel_size*1e6 , '.')
    title('Vertical shift and FOV')
    axis tight
    ylabel('Estimated shift [um]')
    xlabel('Angle [deg]')
    grid on 

    plotting.suptitle('Shifts estimated from initial low-res tomogram')
    
    drawnow 
    
    utils.verbose(struct('prefix', 'template'))

    
end


function [stack_object, total_shift] = find_alignment(stack_object, tomo_lres,theta, par, resolution_ratio, win_small, Nw_lres)
    import utils.* 
    import math.* 
    

    Npix_lres = size(tomo_lres);
    Nw_local = [size(stack_object,1), size(stack_object,2)]; 
    
    [cfg_lres, vectors_lres] = astra.ASTRA_initialize(Npix_lres,Nw_lres,theta, par.lamino_angle, 0, 1);

    % find optimal split of the dataset for given GPU 
    split = astra.ASTRA_find_optimal_split(cfg_lres);
    % forward projection model
    model = astra.Ax_partial(tomo_lres,cfg_lres, vectors_lres,split,'verbose', 0);
    model = exp(1i*model);

    % get phase difference
    diff_model =  math.get_phase_gradient_1D(model, 2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% FIND RELATIVE SHIFT BETWEEN THE DOWNSCALED MODEL AND THE LOW RES RECONSTRUCTION %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    win_small = single(win_small)/255; 
      
    object_small = interpolateFT_centered(stack_object,ceil(Nw_local/resolution_ratio), -1);
    object_small = crop_pad(object_small ./ (abs(object_small) + 1e-3) .* win_small,Nw_lres);

    % get at least some initial phase ramp removal 
    [object_small, gamma_tot, gamma_tot_x, gamma_tot_y] = stabilize_phase(object_small, model, 'weights', abs(object_small));

    % find relative shift of the patch and low resolution model
    total_shift = 0; 

    W = abs(object_small); 
    % perform crosscorrelation between phase derivatives to find
    % the optimal shift 
    shift = find_shift_fast_2D( W.* diff_model, W.* math.get_phase_gradient_1D(object_small, 2));
    % shift to the estimated position 
    object_small = imshift_fft(object_small, shift);
    total_shift = total_shift + shift; 


    %% remove phase ramp using the low res tomogram
    for ii = 1:5
        % refine the phase ramp, we need high precision -> do several iterations 
        [object_small, gamma, gamma_x, gamma_y] = stabilize_phase(object_small, model,'weights', abs(object_small), 'fourier_guess', false);
        gamma_tot = gamma_tot .* gamma; 
        gamma_tot_x = gamma_tot_x + gamma_x; 
        gamma_tot_y = gamma_tot_y + gamma_y; 
    end
    % apply results from low resolution to the full resolution object 
    stack_object = apply_ramp(stack_object,gamma_tot, (gamma_tot_x)/resolution_ratio, (gamma_tot_y)/resolution_ratio );
    total_shift = total_shift * resolution_ratio; 
    
%     weight = uint8(255*win_small);
    
%     weight = uint8(255*real(abs(object_small)));
%     weight(weight<10) = 0;   % remove interpolation artefacts from regions far from measured array 

end


function object_full = apply_ramp(object_full,gamma, gamma_x, gamma_y )
    [M,N,~] = size(object_full);
    xramp = pi*(linspace(-1,1,M))';
    yramp = pi*(linspace(-1,1,N));
    if ~isa(object_full, 'gpuArray')
        object_full = bsxfun(@times,object_full , gamma); 
        object_full = bsxfun(@times,object_full , exp(1i*bsxfun(@times,xramp, M*gamma_x))); 
        object_full = bsxfun(@times,object_full , exp(1i*bsxfun(@times,yramp, N*gamma_y))); 
    else
        % use inplace GPU calculation
        object_full = arrayfun(@auxfun, object_full, gamma, M*gamma_x, N*gamma_y, xramp, yramp);
    end
end

function  object = auxfun(object, gamma, gamma_x, gamma_y, xramp, yramp)

    object = object .* gamma; 
    object = object .* exp(1i*xramp*gamma_x); 
    object = object .* exp(1i*yramp*gamma_y); 

end

