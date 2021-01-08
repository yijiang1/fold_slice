%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example script for ASTRA wrappers 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% recompile commands
%  (Linux, GCC 4.8.5)   mexcuda -outdir private  ASTRA_GPU_wrapper/ASTRA_GPU_wrapper.cu ASTRA_GPU_wrapper/util3d.cu ASTRA_GPU_wrapper/par3d_fp.cu ASTRA_GPU_wrapper/par3d_bp.cu
%  (Windows)  mexcuda -outdir private  ASTRA_GPU_wrapper\ASTRA_GPU_wrapper.cu ASTRA_GPU_wrapper\util3d.cu ASTRA_GPU_wrapper\par3d_fp.cu ASTRA_GPU_wrapper\par3d_bp.cu


% volume settings 
Npix_vol = [300, 300, 100] ;
Nangles = 400; 
Npix_proj = [400, 400]; 

% create "data"
angles = linspace(0, 360, Nangles); 
lamino_angle = 60; 
volData = ones(Npix_vol, 'single'); 

% generate geometry 
[cfg, vectors]  = astra.ASTRA_initialize(Npix_vol, Npix_proj, angles, lamino_angle); 
% find optimal split, for small volumes below 600^3 no split is needed 
split  = astra.ASTRA_find_optimal_split(cfg); 

% generate projections 
projData  = astra.Ax_partial(volData, cfg, vectors, split); 

figure(1)
% plot the projections 
subplot(1,2,1)
plotting.imagesc3D(projData); axis off image; colormap bone
title('Angular geometry')

% do backprojection projections !!! not FBP !!!
backprojData  = astra.Atx_partial(projData, cfg, vectors, split); 


% generate rotation matrix, note that the expected angles needs to be
% adjusted to provide same results are the previous example 
R3 = utils.get_rotation_matrix_3D((90-lamino_angle)*ones(Nangles,1), -angles, zeros(Nangles,1)); 

% generate geometry 
[cfgR, vectorsR]  = astra.ASTRA_initialize(Npix_vol, Npix_proj, R3); 

% generate projections 
projDataR  = astra.Ax_partial(volData, cfgR, vectorsR, split); 

figure(1)
% plot the projections 
subplot(1,2,2)
% plot the projections 
plotting.imagesc3D(projDataR); axis off image; colormap bone
title('Rotation matrix geometry')




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
