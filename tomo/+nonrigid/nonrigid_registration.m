% nonrigid_registration - recovered DVF for given full reconstruction and
% subreconstructions using optical flow method 
%
%[shift_3D_total, vol_err, img_deform] = nonrigid_registration(volume_deform, volume_reference, weight, par, smooth, Niter, shift_3D_total)
%
% Inputs:
%    **volume_deform        low quality deformed volume to be matched with reference 
%    **volume_reference     reference volume used for alignment 
%    **weight               importance weights for the 3D volumes 
%    **par                  tomography parameter structure 
%    **smooth               constant to smooth the recovered DVF 
%    **Niter                number of iterations for the optical flow method 
%    **shift_3D_total       initial deformation field (zeros) 
% Outputs: 
%    ++shift_3D_all         recovered deformation for given full reconstruction and subreconstructions 
%    ++vol_err              error between volume and subvolumes 

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

function [shift_3D_total, vol_err] = nonrigid_registration(volume_deform, volume_reference, weight, par, smooth, Niter, shift_3D_total)
    % estimate deformation fields to match two 3D volumes 

    Npix = size(volume_deform);
        
    if ~exist('shift_3D_total', 'var')
        for ii = 1:3
            shift_3D_total{ii}= gpuArray.zeros( ceil(Npix/par.downsample_DVF) , 'single');
        end
    end
    
    img_deform =  utils.interp3_gpu(volume_deform, shift_3D_total{:});
    
    for iter = 1:Niter

        % core : estimate of the deformation field 
        [shift_3D,vol_err(iter)] = nonrigid.find_shift_3D_nonrigid(img_deform,volume_reference, weight, par.downsample_DVF, smooth, par.regular);
        for ii = 1:3
            shift_3D_total{ii} = shift_3D_total{ii}+ par.relax_pos_corr*shift_3D{ii};
        end
        % apply inverse deformations on the deformated object 
        img_deform =  utils.interp3_gpu(volume_deform, -shift_3D_total{1}, -shift_3D_total{2},  -shift_3D_total{3} );
        
        if iter > 1 && vol_err(end) > vol_err(end-1)
            break
        end
    end

end