%  UNWRAP_2D_BOOTSTRAP Refine sinogram using tomography self-consitency ->
%  try to improve reconstruction if the phase-gradients are too large or
%  dataset contain residua  and other unwrapping methods do not work well. 
%  It is computationally significantly slower than utils.unwrap_2D methods
%  
%  METHOD: 
%  This methods reconstructs tomogram in 2x lower resolution to gain 
%  "redundancy" between the projections. Then synthetic projection of this tomogram 
%  are subtracted from the measured complex projections -> P_difference = P_orig * conj(-i*phase_synthetic_unwrapped)
%  and updated phase is estimated as phase_n = phase_(n-1) + unwrap_2D(P_difference)
%  This bootstrap procedure is repeated in several iteratios. If |P_difference| < pi in some projections 
%  exact unwrapping using phase_n = phase_(n-1) + angle(P_difference) is used. 
%
%  [sinogram] = unwrap_2D_bootstrap(object, theta ,par, Niter) 
%
% Inputs:
%     **object    - complex valued projections 
%     **theta     - initial sinogram guess 
%     **par       - ASTRA config file 
%     **Niter     - ASTRA config vectors 
% Outputs: 
%     ++sinogram  - improved unwrapping of the phase sinogram 

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


function [sinogram] = unwrap_2D_bootstrap(object, theta ,par, Niter, ROI) 
    % try to refine the sinogram using FBP reconstruction as intial guess 
    import utils.* 
    import math.* 
    
    binning = 2; 
    method = 'FBP'; 

    utils.verbose(struct('prefix', 'unwrap'))
   
    % important for laminography case 
%     weights = tomo.Ax_sup_partial(ones([Npix,Npix,Nlayers], 'single'), cfg, vectors,[1,1,Ngpu],tomo_params{:});
%     weights = gather(weights / max(weights(:))); 
%         

    verbose(0,'Bootstrap unwrapping')
    % get initial 2D-FFT phase unwrapping
    sinogram = -tomo.unwrap2D_fft2_split(object,par.air_gap,0,[],par.GPU_list,ROI); 
    sinogram_0 = sinogram; 

    verbose(0,'2D downsampling')
    Np = size(sinogram); 
    sinogram_small = tomo.block_fun(@utils.interpolateFT_centered,sinogram,ceil(Np(1:2)/2/binning)*2, -1); 


        
    [Nlayers,width_sinogram,~] = size(sinogram_small);
    Npix = ceil(width_sinogram/sqrt(2)/32)*32;  % for pillar it can be the same as width_sinogram; 
    [cfg, vectors] = astra.ASTRA_initialize([Npix,Npix, Nlayers],[Nlayers,width_sinogram],theta,par.lamino_angle); 
    % find optimal split of the dataset for given GPU 
    Ngpu = max(1,length(par.GPU_list)); 
    split = astra.ASTRA_find_optimal_split(cfg, Ngpu);
    tomo_params = { 'split', [1,1,Ngpu*split(3)], 'split_sub',[split(1:2),1],  'GPU', par.GPU_list, 'verbose', 1}; 

    residua = tomo.block_fun(@aux_get_residua,object);
    if all(residua == 0)
        verbose(0,'No residua detected, returning FFT_2D unwrapping result')
        [sinogram] = tomo.block_fun(@update_sinogram,object, sinogram_small, par,binning, struct('ROI', {ROI}));
        return
    end
    

    for ii = 1:Niter
        switch method
            case 'CGLS'
                verbose(0,'CGLS')
                rec  = tomo.CGLS(rec, sinogram_small, cfg, vectors, Niter_tomo, tomo_params{:});
            case 'FBP'
              verbose(0,'FBP')
              rec  =  tomo.FBP_zsplit(sinogram_small, cfg, vectors,tomo_params{:});
        end
                 
        % "positivity" constraint
        rec = max(0, rec); 

        verbose(0,'Projection ')
        sinogram_small_updated = tomo.Ax_sup_partial(rec, cfg, vectors, [1,1,Ngpu*split(3)], tomo_params{:});
                
        [sinogram, sinogram_small, upd_norm(ii,:)] = tomo.block_fun(@update_sinogram,object, sinogram_small_updated, par,binning, struct('ROI', {ROI}));
       
        %% plot evolution 
        
        plotting.smart_figure(244)
        subplot(1,2,1)
        plot(mean(upd_norm,2))
        title('Sinogram update norm')
        xlabel('Iteration')
        ylabel('Difference between complex-object and sinogram')
        grid on 
        axis tight 
        subplot(1,2,2)
        [~,ind] = sort(theta); 
        % show only projections with some residuas
        ind = ind(ismember(ind, find(residua))); 
        plotting.imagesc3D(cat(2, sinogram_0(:,:,ind), sinogram(:,:,ind))); 
        title('Original sinogram (left)     Improved sinogram (right)')
        axis off xy image 
        colormap bone 
        plotting.suptitle('Bootstrap unwrapping')
        win_size = [1400 500]; 
        screensize = get( groot, 'Screensize' );
        set(gcf,'Outerposition',[150 min(270,screensize(4)-win_size(2)) win_size]);

        drawnow 

    end
    
    utils.verbose(struct('prefix', 'template'))


end

function [sinogram, sinogram_small, upd_norm] = update_sinogram(object, sinogram_small, par, binning)
    Np = size(object); 
    % upsample small sinogram back to the full size
    sinogram = utils.interpolateFT_centered(sinogram_small,Np(1:2), -1); 
   
    % use the knowledge that around phase jumps is usually zero or very low intensity
    W = min(1, abs(object)); 
    
    
    %% sinogram refinement 
    % find sinogram ramp and offset to match the tomo guess
    object_resid = object.*exp(1i*sinogram); 

    % estimate the update using 2D phase unwrap 
    sinogram = sinogram - W.*math.unwrap2D_fft2(object_resid,par.air_gap,0);
 
    % make sinogram exactly equal to the data , 
    % !! dangerous, it can make it even worse 
    % -> allow it only for the well behaved projections 
    phase_update = angle(object.*exp(1i*sinogram)); 
    minor_update_ind = all(all(abs(phase_update)<0.5)); 
    sinogram = sinogram - minor_update_ind.*W.*phase_update; 


    upd_norm = squeeze(math.norm2(angle(object_resid)));
    
    % get a downsampled version of the sinogram 
    sinogram_small = utils.interpolateFT_centered(sinogram,ceil(Np(1:2)/2/binning)*2, -1); 

end

function residua = aux_get_residua(object_block)
    % GPU auxiliarly function 
    residua = squeeze(math.sum2(abs(utils.findresidues(object_block))>0.1));
end
