% FBP_PROPAGATION filtered back propagation for diffraction tomography 
%
% [rec] = FBP_propagation(sino, theta, variable, par, optimal_propagation)
%
% Inputs:
%     **sino        - sinogram (Nlayers x width x Nangles)
%     **angles      - projection angles 
%     **variable    - 'phase' or 'amplitude'
%     **par         - parameter structure 
%     **optimal_propagation - position of center of focus 
% *returns*
%      ++rec - reconstructed volume 

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

function [rec_volume] = FBP_propagation(sino, theta, variable, par, optimal_propagation, thickness)
% Filtered backpropagation

% assert(~isreal(sino), 'Input sinogram has to be complex-valued')

[Nlayers, width_sinogram, Nangles] = size(sino); 

assert(length(theta)==Nangles, 'Size of sinogram does not correspond to size of Nangles');

% create a matrix of propagation through entire sample 
if nargin < 6 
    propag = -single(par.pixel_size*(-ceil(width_sinogram/2):floor(width_sinogram/2-1)));
else
    propag = -single(linspace(-thickness/2,thickness/2, width_sinogram));
end
[~,H0] = utils.prop_free_nf(ones(Nlayers,width_sinogram,'single'), par.lambda, propag, par.pixel_size); 

%% allocate shared memory 
use_sharemem = false; 
if use_sharemem
    %% allocate reconstruction volume 
    rec_volume = zeros(width_sinogram, width_sinogram, Nlayers, 'single'); 
    share_mem = shm(true);
    share_mem.allocate(rec_volume); 
    share_mem.detach(); 
else
    %% allocate reconstruction volume 
    rec_volume = gpuArray.zeros(width_sinogram, width_sinogram, Nlayers, 'single'); 
    H0 = gpuArray(H0);
end

% create a support mask that limits extend of the reconstruction 
[~,circle] = utils.apply_3D_apodization(rec_volume,50,0,0.1); 
circle = single(circle); 


%% solve the FBP task 
for ii = 1:Nangles
    utils.progressbar(ii, Nangles, 100); 
    rec_volume = calculate_filt_back_propagation(rec_volume, sino(:,:,ii), theta(ii),circle,H0, par, variable, optimal_propagation(min(ii,end))); 
end

if use_sharemem
    [share_mem, rec_shm] = share_mem.attach(); 
    rec_volume(:) = rec_shm; 
    share_mem.detach(); 
end
    
rec_volume = gather(rec_volume); 

end

function rec_volume = calculate_filt_back_propagation(rec_volume, sino_tmp, theta, circle, H0, par, variable, optimal_propagation, thickness)
    
    [Nlayers, width_sinogram] = size(sino_tmp); 
    
    sino_tmp = propagate_sinogram(H0, sino_tmp, par,variable, optimal_propagation);
    
    sino_tmp = permute(sino_tmp, [3,2,1]);

    cfg.iProjV = size(sino_tmp,1);
    cfg.iProjU = width_sinogram;
    cfg.iProjAngles = Nlayers;

    % apply filtering, dont do any backprojetion 
    [~,sino_filt] = tomo.FBP(sino_tmp, cfg, zeros(Nlayers,12), 1,'verbose',0, 'determine_weights', false, 'GPU', par.GPU_list, 'filter', 'ram-lak', 'only_filter_sinogram', true);

    
    if isscalar(H0)
        sino_filt = repmat(sino_filt,width_sinogram,1,1); 
    end
        
    sino_filt = sino_filt .* circle;

    %% rotate propagated projections 
    % sino_filt = utils.imrotate_ax_fft(sino_filt, theta(ii), 3); % subpixel precision inteprolation using FFT
    sino_filt = utils.imrotate_ax(sino_filt, theta, 3);  % common linear interpolation 
    
    %% add filtered update to the total reconstruction 
    if isa(rec_volume, 'shm')
        [share_mem, rec_volume] = rec_volume.attach(); 
        tomo.set_to_array(rec_volume, gather(sino_filt), 0, true);  % write directly to the shared memory 
        share_mem.detach();
    else
        rec_volume = rec_volume + sino_filt; 
    end
    

end
        
function sinogram = propagate_sinogram(H0, sinogram, par,variable, optimal_propagation)

    [Nlayers, width_sinogram, ~] = size(sinogram); 
    
    sinogram = gpuArray(sinogram); 
    
    
    
    %% FT interpolation + NF propagation 
    if any(optimal_propagation > 0)
        H = gpuArray.ones(Nlayers,width_sinogram,'single'); 
        [~,H] = utils.prop_free_nf(H, par.lambda, optimal_propagation, par.pixel_size); 
    else
        H = 1; 
    end
    if ~isscalar(H0) || ~isscalar(H)
        % propagate along the beam 
        sinogram = ifft2(fft2(sinogram).*H0.* H); 
    end
    
    switch variable 
        case 'amplitude' 
        %% get amplitude 
            sinogram = -log(abs(sinogram)); 
        case 'phase'
        %% get phase
            sinogram = -math.unwrap2D_fft(sinogram,2, [10,10], 0); 
        otherwise
            error('Wrong option')
    end
    
end
