function [dof,pixel_size] = DOF_calculator(energy, det_pixel, det_N, distance, alpha)
%Calculate theoreical depth of focus for X-ray ptychography

% Inputs:  
%     **energy     beam energy (keV)
%     **det_pixel  detector pixel size (m) 
%     **det_N      # of pixels in the detector
%     **distance   # sample to detector distance (m)
%     **alpha      additional scaling coefficient
% *returns*: 
%     ++dof    depth of focus for ptychography

lambda = 1.23984193e-9/energy; % wavelength (m)
pixel_size = lambda*distance/(det_pixel)/det_N; %pixel size in ptycho reconstruction (m)
dof = alpha * pixel_size^2/lambda;
end

