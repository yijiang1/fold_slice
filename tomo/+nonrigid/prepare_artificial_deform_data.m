% prepare_artificial_deform_data - prepare deformed phantom for algorithm tests
%
%[dphase, shift_3D_orig, angles, par, rec_ideal, volData_orig] = ...
%    prepare_artificial_deform_data(Nangles, Npix, Nlayers, Nblocks, smooth, binning,DVF_amplitude,DVF_period,  par_0)
%
% Inputs:
%    **Nangles              number of angles in the simulated dataset 
%    **Npix                 int - pixel size of the phantom 
%    **Nlayers              number of layers in the phatom 
%    **Nblocks              number of subtomograms 
%    **smooth               constant used to estimate ratio between phantom pixel size and DVF pixels size            
%    **binning              simulate binning of the produced sinograms 
%    **par_0                initial parameters structure that will be merged with the loaded paramters 
% Outputs: 
%    ++dphase               phase difference for the complex project 
%    ++shift_3D_all_0       original deformation = {}
%    ++angles               angles for each of the projection 
%    ++par                  merged parameter structure 
%    ++rec_ideal            ideal construction with known DVF
%    ++volData_orig         original phantom 

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


function [dphase, shift_3D_orig, angles, par, rec_ideal, volData_orig] = ...
    prepare_artificial_deform_data(Nangles, Npix, Nlayers, Nblocks, smooth, binning,DVF_amplitude,DVF_period, par_0)

%% prepare deformated data for provided parameters 

import utils.*
import math.*
import plotting.*

%% create data and geometry
angles = pi+[linspace(0, 180, Nangles)];
% create 8 subtomos 
angles = reshape(angles, Nblocks,[])';
angles = angles(:); 

lamino_angle = 90; 
Nw = ceil(Npix*sqrt(2)/16)*16;


try
    disp('Loading stored model')
    load(['porous_glass_data',num2str(Npix),'.mat']); 
    disp('Loading done')
catch
    try
        %% porous_glass phantom 
        disp('Creating phantom')
        rng default 
        volData_orig = randn(2*[Npix, Npix, Nlayers], 'single');
        volData_orig = imgaussfilt3_fft(volData_orig, 6);
        porous_glass = imgaussfilt3_fft(volData_orig > 0, 2)>0.01 & volData_orig <= 0;
        volData_orig = randn(2*[Npix, Npix, Nlayers], 'single');
        volData_orig = imgaussfilt3_fft(volData_orig, 6);
        porous_glass = porous_glass | imgaussfilt3_fft(volData_orig > 0, 2)>0.01 & volData_orig <= 0;
        [Xq,Yq,Zq] = meshgrid(linspace(-0.5,0.5,2*Npix), linspace(-0.5,0.5,2*Npix), linspace(-0.5,0.5,2*Nlayers));
        porous_glass = interp3(single(porous_glass),Xq*Npix*3.5+Npix,Yq*Npix*3.5+Npix,Zq*Nlayers*3.5+Nlayers); 
        porous_glass = porous_glass(end/4:end*3/4-1,end/4:end*3/4-1,end/4:end*3/4-1);
        porous_glass(isnan(porous_glass)) = 0; 

        % apply circular mask 
        xgrid = -Npix/2+1 : Npix/2;
        [X,Y]  = meshgrid(xgrid, xgrid);
        porous_glass = porous_glass .* imgaussfilt(single(X.^2+Y.^2 < (Npix/2.2)^2), 3);
        porous_glass = porous_glass .* reshape(tukeywin(Nlayers, 0.5), 1,1,[]);
        porous_glass = uint8(porous_glass/max(porous_glass(:)) * 255); 
        savefast_safe(['porous_glass_data',num2str(Npix),'.mat'], 'porous_glass', true);
    catch
        keyboard
    end
end



Bsize = ceil(Nangles/Nblocks); 

% load porous_glass_data
volData_orig = single(porous_glass); 
volData_orig = volData_orig(:,:,1:Nlayers);

% apply circular mask 
xgrid = -Npix/2+1 : Npix/2;
[X,Y]  = meshgrid(xgrid, xgrid);
volData_orig = volData_orig .* imgaussfilt2_fft(single(X.^2+Y.^2 < (Npix/2.4)^2), 5);
volData_orig = volData_orig .* reshape(tukeywin(Nlayers, 0.1), 1,1,[]);
Nlayers = size(volData_orig,3);

%% initialize deformation vector fields reconstructions 
Nps = ceil([Npix, Npix, Nlayers]/par_0.downsample_DVF);


%% generate deformation field 
volData_orig = gather(volData_orig);


%% generate "measured" data 

disp('Generating data')
split = 1;

rng default 
for ax= 1:3
    for j = 1:2
        shift_3D{j}{ax} = imgaussfilt3_fft(randn(Nps), DVF_period);
        shift_3D{j}{ax} = shift_3D{j}{ax} / max(abs(shift_3D{j}{ax}(:)))*DVF_amplitude;
    end
end

for ll = 1:Nblocks+1
%     ratio(1) =  1-exp(-3*((ll-1)/(Nblocks+1)));
%     ratio(2) =  1-exp(-3*((ll-1)/(Nblocks+1)));
%     ratio(3) =  1-exp(-3*((ll-1)/(Nblocks+1)));
    
    ratio(1) =  sin(2*pi*(ll-1)/(Nblocks+1));
    ratio(2) =  sin(2*pi*(ll-1)/(Nblocks+1));
    ratio(3) =  sin(2*pi*(ll-1)/(Nblocks+1));
        
    for ax= 1:3
        shift_3D_orig{ll}{ax} = (ratio(ax)*shift_3D{1}{ax});
    end 
end


[cfg, vectors] = ...
astra.ASTRA_initialize([Npix, Npix,Nlayers],[Nlayers,Nw],angles,lamino_angle, 0, 1);


% resample the created DVF to reconstruction size of the DVF 
for ll = 1:Nblocks+1
    for kk  = 1:3
        Np = size(shift_3D_orig{ll}{kk}); 
        [X,Y,Z] = meshgrid(linspace(1,Np(1),Nps(1)), linspace(1,Np(2),Nps(2)), linspace(1,Np(3),Nps(3)));
        shift_3D_orig{ll}{kk} = interp3(shift_3D_orig{ll}{kk},X,Y,Z);
    end
end

[deform_tensors,inv_deform_tensors] = nonrigid.invert_DVF(shift_3D_orig, [Npix, Npix,Nlayers]);
% join blocks to keep initial and final deform for each block together 
for block = 1:Nblocks
    deform_tensors_linear{block} = [deform_tensors{block}; deform_tensors{block+1}]; 
    inv_deform_tensors_linear{block} = [inv_deform_tensors{block}; inv_deform_tensors{block+1}]; 
end



sinogram = tomo.Ax_sup_partial(volData_orig,cfg, vectors,split);

rec_ideal = tomo.FBP(sinogram , cfg, vectors,split, 'verbose',0);   

% generate data 
for ll = 1:Nblocks
    ids = 1+(ll-1)*Bsize:min(Nangles, ll*Bsize);
    cfg.iProjAngles = length(ids);
    
    sinogram(:,:,ids) = tomo.Ax_sup_partial(volData_orig,cfg, vectors(ids,:),split, ...
        'deformation_fields', deform_tensors_linear{ll});
 
end


%create realistic issues 
sinogram = binning_2D(sinogram,binning); 

% change the change to get phase jumps 
sinogram = sinogram / max(sinogram(:)) * 2*pi;

dphase = math.get_phase_gradient_1D(-sinogram, 2);



[cfg, vectors] = ...
astra.ASTRA_initialize([Npix, Npix,Nlayers]/binning,[Nlayers,Nw]/binning,angles,lamino_angle, 0, 1);


rec_0  = tomo.FBP(sinogram , cfg, vectors,split);
rec_corr = 0; 
for ll = 1:Nblocks
    ids = 1+(ll-1)*Bsize:min(Nangles, ll*Bsize);
    cfg.iProjAngles = length(ids);
    rec_corr = rec_corr+tomo.FBP(sinogram(:,:,ids) , cfg, vectors(ids,:),split, 'verbose',0,...
        'deformation_fields', inv_deform_tensors_linear{ll} )/Nblocks;
end


% if debug()
    figure
    subplot(1,2,1)
    imagesc3D(max(0,rec_0), 'init_frame', Nlayers/2)
    axis off image; colormap bone 
    title('Standard reconstruction')
    subplot(1,2,2)
    imagesc3D(max(0,rec_corr), 'init_frame', Nlayers/2)
    axis off image; colormap bone 
    title('Ideally corrected reconstruction')
    drawnow 
% end

% store inputs to par structure 
par.binning = binning; 
par.valid_angles = 1:Nangles;
par.air_gap = [20,20];
par.factor = 1; 
par.output_folder = ''; 

for field = fields(par_0)'
    par.(field{1}) = par_0.(field{1}); 
end





end


