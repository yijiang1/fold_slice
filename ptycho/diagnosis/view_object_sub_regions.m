close all
addpath(fullfile(pwd,'utils'))

%% load a reconstruction
disp('extracting sub-regions...')
[oROI, oROI_vec, sub_px_shift] = find_reconstruction_ROI_external( outputs.probe_positions, size(object), size(probe(:,:,1,1)) );
object_sub_regions = zeros(size(probe,1),size(probe,2),size(outputs.probe_positions,1));
for i=1:size(object_sub_regions,3)
    object_sub_regions(:,:,i) = object(oROI_vec{i,:});
end
disp('extracting sub-regions...done')

%% show in sub regions of object phase
imagesc3D(angle(object_sub_regions))
colormap gray
axis image

%% save sub regions
disp('saving sub regions...')
saveName = 'Niter1000_sub_regions_ph.hdf5';
N = size(object_sub_regions,1);
h5create(saveName, '/sub_recons_phase', size(object_sub_regions),'ChunkSize',[N N, 1],'Deflate',4)
h5write(saveName, '/sub_recons_phase', angle(object_sub_regions))
saveName = 'Niter1000_sub_regions_mag.hdf5';
h5create(saveName, '/sub_recons_mag', size(object_sub_regions),'ChunkSize',[N N, 1],'Deflate',4)
h5write(saveName, '/sub_recons_mag', abs(object_sub_regions))
disp('saving sub regions...done')

%%
saveName = 'Niter1000_positions.hdf5';
ppX = outputs.probe_positions(:,1)*p.dx_spec(1);
ppY = outputs.probe_positions(:,2)*p.dx_spec(1);
hdf5write(saveName, '/ppX', ppX)
hdf5write(saveName, '/ppY', ppY,'WriteMode','append')
hdf5write(saveName, '/dx', p.dx_spec(1),'WriteMode','append')
