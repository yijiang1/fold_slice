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
figure
imagesc3D(abs(psi_s))
colormap parula
axis image

%% Generate diffraction patterns from reconstructed object and probes
% Note: this is a very crude approximation because
% 1. The layer distance could be too large for accurate forward propagations.
% 2. Only the 1st probe mode is used.
% 3. Sub-pixel shifts are NOT considered.

Np_p = size(probe(:,:,1,1));
N_scans = size(outputs.probe_positions,1);
Nlayer = length(p.multi_slice_param.z_distance)-1;
dp_fwd = zeros(Np_p(1), Np_p(1), N_scans);

% calculate propagator
[~,H,~,~] = near_field_evolution(ones(Np_p), p.multi_slice_param.z_distance(1), p.lambda,  p.dx_spec.*Np_p, true );
H = ifftshift(H);

% get object indicies for each scan point
[oROI, oROI_vec, sub_px_shift] = find_reconstruction_ROI_external(outputs.probe_positions, size(object(:,:,1)), size(probe(:,:,1,1)) );
object_sub_regions = zeros(size(probe,1), size(probe,2), size(outputs.probe_positions,1));

% if OPR is used, calculate probes at each scan position
if size(probe,4)>1
    probes = reshape(probe(:,:,1,:), prod(Np_p), []);
    probes = reshape(probes * outputs.probe_evolution', Np_p(1), Np_p(2), []);
end

% calculate diffraction patterns at each scan position
f = waitbar(0, '1', 'Name', 'Generating diffraction patterns...',...
        'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
setappdata(f, 'canceling', 0);

for i=1:N_scans
    % Check for clicked Cancel button
    if getappdata(f, 'canceling')
        break
    end
    % Update waitbar and message
    waitbar(i/N_scans, f, sprintf('No.%d/%d', i, N_scans))

    if size(probe, 4)>1
        psi = probes(:,:,i);
    else
        psi = probe(:,:,1,1);
    end
    
    for k=1:Nlayer
        psi = psi.*object(oROI_vec{i,:}, k);
        psi = ifft2(bsxfun(@times, H, fft2(psi)));
    end
    dp_fwd(:,:,i) = abs(fftshift(fft2(psi))).^2;
end
delete(f)

figure
imagesc3D(dp_fwd)
colormap jet
axis image

%% compare raw data and simulated diffraction patterns
% load raw data
dp_raw = h5read('data_roi0_Ndp256_dp.hdf5','/dp');
% additional processing (e.g., crop, flip) might be needed

figure
dp_combined = cat(2, dp_raw, dp_fwd);
imagesc3D(log(dp_combined))
%imagesc3D(dp_combined)

colormap jet
axis image
