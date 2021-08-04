%view_multi_slice_recon.m
addpath(fullfile(pwd,'utils'))

%load a matlab recon file

%% View object phase
object_roi = object(p.object_ROI{:},:);
object_roi_unwrapped = object_roi;
for i=1:size(object_roi,3)
    object_roi_unwrapped(:,:,i) = phase_unwrap(angle(object_roi(:,:,i)));
end

imagesc3D(object_roi_unwrapped)
colormap gray
axis image

%% View probe propagation
% calculate propagator
Np_p = size(probe(:,:,1,1));
[~,H,~,~] = near_field_evolution(ones(Np_p), p.multi_slice_param.z_distance(1), p.lambda,  p.dx_spec.*Np_p, true );
H = ifftshift(H);
Nlayer = length(p.multi_slice_param.z_distance)-1;
psi_s = zeros(Np_p(1),Np_p(2),Nlayer);
psi_s(:,:,1) = probe(:,:,1,1);
for k=2:Nlayer
    psi_s(:,:,k) = ifft2(bsxfun(@times, H, fft2(psi_s(:,:,k-1))));
end