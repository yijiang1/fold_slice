function [mask] = get_tomo_fourier_mask_3d( N_recon, angles )

mask = tomo.get_tomo_fourier_mask_2d(N_recon(2),[N_recon(1),N_recon(2)],angles);
mask = repmat(mask,[1,1,N_recon(3)]);
%mask = mask==0;
mask = ifftshift(mask);
end