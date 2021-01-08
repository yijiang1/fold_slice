function [recon_constrained] = TV_Fourier_smoothing(I, I_f, Niter, show_figure)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%Convergence Parameters
%Niter = 10;
iter_TVdecent = 30; 
a = 0.5; %Decent Parameter

% Create a mask that will remove the region of interest (ROI). 
mask = I_f~=0;

%Create Random Image
[nx, ny] = size(I);
recon_init = rand(nx,ny);

%recon_init = I;

%% 
for i = 1:Niter
    
    % Counter. 
    disp(i)
          
    % FFT of Reconstructed Image.
    FFTr = fftshift(fft2(ifftshift(recon_init))); 
    
    % Remove the ROI with Data Constraint. 
    FFTr(mask) = I_f(mask);
 
    %Inverse FFT
    recon_constrained = real(fftshift(ifft2(ifftshift(FFTr))));
    
    %Positivity Constraint.
    %recon_constrained(recon_constrained<0) = 0;

    %TV Minimization. 
    recon_minTV = recon_constrained;
    d = (sum(sum((recon_minTV-recon_init).^2))).^(1/2);
    for j = 1:iter_TVdecent

        Vst = TVDerivative(recon_minTV);
        L2norm = (sum(sum(Vst.^2))).^(1/2); 
        Vst = Vst/L2norm;
        recon_minTV = recon_minTV - a*d*Vst;

    end

    % Initialize next loop. 
    recon_init = recon_minTV;
    
end
if show_figure
    % Show the Reconstruction
    figure
    imagesc((recon_minTV+recon_constrained)/2); axis image; colormap gray
    figure
    f = fftshift(fft2(recon_constrained));
    imagesc(abs(f).^0.2); axis image; colormap jet
end
%Save the Reconstructions
%imwrite(mat2gray(recon_constrained), [fname '_Reconstruction.tif'], 'tiff') 

end

