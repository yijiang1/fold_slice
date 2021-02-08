%view_multi_slice_recon.m
addpath(fullfile(pwd,'utils'))

%load a matlab recon file
%% unwrap phase
object_roi = object(p.object_ROI{:},:);
object_roi_unwrapped = object_roi;
for i=1:size(object_roi,3)
    object_roi_unwrapped(:,:,i) = phase_unwrap(angle(object_roi(:,:,i)));
end
%%
imagesc3D(object_roi_unwrapped)
colormap gray
axis image