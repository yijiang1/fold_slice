% show_reconstruction_quality - compare conventional reconstruction,
% nonrigid reconstruction FBP and SART 
%
%   [rec_FBP, rec_NCT_FBP, rec_NCT_SART] = show_reconstruction_quality(sinogram, cfg, vectors,  shift_3D_total, regularize_deform_evol)
%
% Inputs:
%    **sinogram          current reconstruction 
%    **vectors          ASTRA configuration vectors
%    **cfg              ASTRA configuration structure
%    **shift_3D_total   recovered deformation vector field 
%    **regularize_deform_evol    regularization constant for the deformation field evolution calculation 
%
% Outputs: 
%    ++rec_FBP          conventional FBP reconstruction
%    ++rec_NCT_FBP      nonrigid FBP reconstruction
%    ++rec_NCT_SART     nonrigid SART reconstruction
% 


%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2018 by Paul Scherrer Institute (http://www.psi.ch)    |
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



function [rec_FBP, rec_NCT_FBP, rec_NCT_SART] = show_reconstruction_quality(sinogram, cfg, vectors,  shift_3D_total, regularize_deform_evol)

Nangles = cfg.iProjAngles; 

reset(gpuDevice)

Nblocks = length(shift_3D_total); 
split = astra.ASTRA_find_optimal_split(cfg);
Bsize = ceil(Nangles/Nblocks); 

if any(split(1:3)) > 1
    warning('Sample volume seems too large, try to reduce the reconstructed volume size')
    split(1:3) = 1;  % at least try to make it work without splitting, otherwise recontruction will be poor 
end

rec_FBP = gather(tomo.FBP(sinogram, cfg, vectors)) ;

% generate deformation tensors from the shift tensors 
[deform_tensors, inv_deform_tensors] = nonrigid.get_deformation_fields(shift_3D_total, regularize_deform_evol, size(rec_FBP));


%% FBP 
rec_NCT_FBP = nonrigid.FBP_deform(sinogram, cfg, vectors,Bsize, inv_deform_tensors);

    
[SART_cache, cfg_SART] = tomo.SART_prepare(cfg, vectors, Bsize, split);
SART_cache.R = min(1,SART_cache.R);



%% SART - solve it using all constraints
[~,rec_mask] = utils.apply_3D_apodization(rec_NCT_FBP,0); 
rec_NCT_SART = rec_NCT_FBP; 
Niter_SART = 10;
clear err_sart
disp('====== SART ==========')

for kk = 1:Niter_SART
    utils.progressbar(kk, Niter_SART)
    [rec_NCT_SART,err_sart(kk,:)] = tomo.SART(rec_NCT_SART, sinogram, cfg_SART, vectors, SART_cache, split, ...
        'relax',0, 'deformation_fields', deform_tensors,'inv_deformation_fields', inv_deform_tensors, ...
          'constraint', @(x)(max(0,x.*rec_mask)), 'verbose',0); 

%     figure(1343)
%     subplot(1,2,1)
%     plot(err_sart)
%     hold all 
%     plot(mean(err_sart'),'k', 'LineWidth',2)
%     hold off 
%     set(gca, 'xscale', 'log')
%     set(gca, 'yscale', 'log')
%     grid on 
%     axis tight 
%     title('SART error evolution')
%     subplot(1,2,2)
%     plotting.imagesc3D(rec_NCT_SART, 'init_frame', floor(size(rec_NCT_SART,3)/2))
%     axis image 
%     colormap bone 
%     axis off image 
%     drawnow 
    
    
end

% remove edges 
rec_FBP = utils.apply_3D_apodization(rec_FBP, 0); 
rec_NCT_FBP = utils.apply_3D_apodization(rec_NCT_FBP,0); 


figure(10)
if exist('orig_phantom', 'var')  && ~isempty(orig_phantom)
    orig_phantom = utils.crop_pad(orig_phantom, [cfg.iVolX,cfg.iVolY]);
    orig_phantom = orig_phantom ./ mean(orig_phantom(:)) * mean(rec_FBP(:))*0.9; 
    rec_all = gather(cat(2, orig_phantom,rec_FBP, rec_NCT_FBP, rec_NCT_SART)); 
else
    rec_all = gather(cat(2, rec_FBP, rec_NCT_FBP, rec_NCT_SART)); 
end
range = quantile(rec_all(:), [1e-2, 1-1e-2]); 

plotting.imagesc3D(rec_all, 'init_frame', size(rec_all,3)/2)
caxis(range);
axis off image; colormap bone 
title('Original reconstruction / Deform FBP / Deform SART')



end