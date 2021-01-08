% TESTING SCRIPT for ASTRA wrappers 
% This script will run some basic features in the ASTRA wrapper code
% ie reconstruciton , splitting on GPU, splitting between multiple workers 
% compare results with expected values to detect inconsistencies 
% recompile commands
%  (Linux, GCC 4.8.5)   mexcuda -outdir private  ASTRA_GPU_wrapper/ASTRA_GPU_wrapper.cu ASTRA_GPU_wrapper/util3d.cu ASTRA_GPU_wrapper/par3d_fp.cu ASTRA_GPU_wrapper/par3d_bp.cu
%  (Windows)  mexcuda -outdir private  ASTRA_GPU_wrapper\ASTRA_GPU_wrapper.cu ASTRA_GPU_wrapper\util3d.cu ASTRA_GPU_wrapper\par3d_fp.cu ASTRA_GPU_wrapper\par3d_bp.cu

cd(fullfile( fileparts(mfilename('fullpath')), '..'))
addpath('tests')
addpath('utils')
addpath('./')
addpath(find_base_package)
utils.verbose(0)

utils.verbose(0,'Creating data for ASTRA wrapper tests')
if ~exist('GPU_id', 'var'); GPU_id = [1]; end
gpuDevice(GPU_id);

% volume settings 
Npix_vol = [300, 300, 200] ;
Nangles = 400; 
Npix_proj = [400, 400]; 

% create "data"
angles = linspace(0, 180, Nangles); 
lamino_angle = 60;   % 90 deg is normal tomo
tilt_angle = 10;      % rotation in plane of the projection 
CoR_offset = [20,10];   % offset of the center of rotation 
pixel_scale = [1,1] ;     % relative scale of the pixels 
volData = ones(Npix_vol, 'single'); 

% generate geometry 
[cfg, vectors]  = astra.ASTRA_initialize(Npix_vol, Npix_proj, angles, lamino_angle, tilt_angle, pixel_scale,Npix_proj/2+CoR_offset); 
% find optimal split, for small volumes below 600^3 no split is needed 
split  = astra.ASTRA_find_optimal_split(cfg); 

% generate projections 
projData  = astra.Ax_partial(volData, cfg, vectors, split); 

% plot the projections 
figure
plotting.imagesc3D(projData, 'init_frame', 50); axis off image; colormap bone
title('Volume projection')
drawnow 

% do backprojection projections !!! not FBP !!!
backprojData  = astra.Atx_partial(projData, cfg, vectors, split); 

utils.verbose(0,'Simple ASTRA wrapper tested')
%% =============== tests simple "on GPU" splitting ============ 
split = [2,2,2,2]; 

% generate projections 
projData_split  = astra.Ax_partial(volData, cfg, vectors, split); 

err = norm(projData(:) - projData_split(:)) / norm(projData(:)); 
% plotting.imagesc3D(projData - projData_split); colorbar;  axis off image; colormap bone
assert(err < 1e-5, "Splitted projections solver is not fully consistent with the unsplitted one")


% do backprojection projections !!! not FBP !!!
backprojData_split  = astra.Atx_partial(projData, cfg, vectors, split); 
err = norm(backprojData(:) - backprojData_split(:)) / norm(backprojData(:)); 
assert(err < 1e-6, "Splitted backprojections solver is not fully consistent with the unsplitted one")


utils.verbose(0,'On GPU splitting tested')
%% =============== tests additional splitting for tomography ============ 

split = [2,2,2,1]; 
projData_split = tomo.Ax_sup_partial(volData, cfg, vectors, split); 
err = norm(projData(:) - projData_split(:)) / norm(projData(:)); 
if err > 1e-5
    figure
    plotting.imagesc3D(projData - projData_split); colorbar;  axis off image; colormap bone
    title('Projection data difference')
    drawnow
end
if err > 1e-3
    error("Sup-splitted projections solver is not fully consistent with the unsplitted one")
elseif err > 1e-5
    warning("Sup-splitted projections solver is not fully consistent with the unsplitted one, most likely only subpixel errors are present")
end



backprojData_split  = tomo.Atx_sup_partial(projData, cfg, vectors, split); 
err = norm(backprojData(:) - backprojData_split(:)) / norm(backprojData(:)); 
if err > 1e-5
    figure
    plotting.imagesc_tomo(backprojData - backprojData_split); colorbar;  axis off image; colormap bone
    title('Backprojection data difference')
    drawnow 
end
if err > 1e-3
    error("Sup-splitted backprojections solver is not fully consistent with the unsplitted one")
elseif err > 1e-5
    warning("Sup-splitted backprojections solver is not fully consistent with the unsplitted one, most likely only subpixel errors are present")
end

utils.verbose(0,'Extra splitting tested')

%% =============== tests multiGPU solvers  ============ 

utils.verbose(0,'Testing multiGPU solver')

c = parcluster('local');
if exist('local_cluster_jobs', 'dir')
    rmdir('local_cluster_jobs', 's')  % delete folder with jobs (prevent accumulation ) 
end
mkdir('local_cluster_jobs')   % recreate the folder 
c.JobStorageLocation = ['local_cluster_jobs'];


GPU = [1:gpuDeviceCount]; 
split = [2,2,2,2]; 
projData_split = tomo.Ax_sup_partial(volData, cfg, vectors, split, 'GPU', GPU); 
err = norm(projData(:) - projData_split(:)) / norm(projData(:)); 
% plotting.imagesc3D(projData - projData_split); colorbar;  axis off image; colormap bone
assert(err < 1e-5, "Sup-splitted projections solver is not fully consistent with the unsplitted one")


backprojData_split  = tomo.Atx_sup_partial(projData, cfg, vectors, split, 'GPU', GPU); 
err = norm(backprojData(:) - backprojData_split(:)) / norm(backprojData(:)); 
assert(err < 1e-6, "Sup-splitted backprojections solver is not fully consistent with the unsplitted one")

utils.verbose(0,'Multi GPU solver tested')
rmdir('local_cluster_jobs', 's')  % delete folder with jobs (prevent accumulation ) 

utils.verbose(0,'==== All ASTRA wrapper tests passed ===== ')

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



