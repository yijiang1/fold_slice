%view_grid_artifact_removal.m

%% Show phase image and its Fourier magnitude to inspect the effect of grid artifact removal
close all
object_roi_ph = angle(object(p.object_ROI{:}));
%object_roi_ph = angle(object_roi); %old output file

figure    
imagesc(object_roi_ph); axis image; colormap gray
figure
imagesc(abs(fftshift(fft2(object_roi_ph))).^0.2); axis image;

%% Quick way to see the scan step size
%load scan positions
ppY = h5read('data_roi0_Ndp128_para.hdf5','/ppY');
ppX = h5read('data_roi0_Ndp128_para.hdf5','/ppX');
figure
subplot(1,2,1)
plot(diff(ppY),'.')
title('Vertical')
subplot(1,2,2)
plot(diff(ppX),'.')
title('Horizontal')
