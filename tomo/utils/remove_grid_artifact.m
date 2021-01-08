function [I_phase_new, f] = remove_grid_artifact(I_phase, dx, stepSize_x,stepSize_y, windowSize, direction, showFigure)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Nx = size(I_phase,2);
Ny = size(I_phase,1);
dk_x = 1/dx/Nx;
dk_y = 1/dx/Ny;
cen_x = floor(Nx/2)+1;
cen_y = floor(Ny/2)+1;

k_max = 1/dx;
f0 = fftshift(fft2(ifftshift(I_phase)));
f = f0;
dk_s_x = 1/stepSize_x;
dk_s_y = 1/stepSize_y;

switch direction
    case 'xy'
        x_range = ceil(-k_max/2/dk_s_x):floor(k_max/2/dk_s_x);
        y_range = ceil(-k_max/2/dk_s_y):floor(k_max/2/dk_s_y);
    case 'x'
        x_range = ceil(-k_max/2/dk_s_x):floor(k_max/2/dk_s_x);
        y_range = 0;
    case 'y'
        x_range = 0;
        y_range = ceil(-k_max/2/dk_s_y):floor(k_max/2/dk_s_y);
end
for i=1:length(x_range)
    for j=1:length(y_range)
        
        if ~(x_range(i)==0 && y_range(j)==0)
            window_x_lb = max(round(x_range(i)*dk_s_x/dk_x) + cen_x - windowSize, 1);
            window_x_ub = min(round(x_range(i)*dk_s_x/dk_x) + cen_x + windowSize, Nx);
            
            window_y_lb = max(round(y_range(j)*dk_s_y/dk_y) + cen_y - windowSize, 1);
            window_y_ub = min(round(y_range(j)*dk_s_y/dk_y) + cen_y + windowSize, Ny);
            %window_y_lb
            %window_y_ub
            f(window_y_lb:window_y_ub,window_x_lb:window_x_ub) = 0;
        end
    end
end
I_phase_new = real(fftshift(ifft2(ifftshift(f))));

if showFigure 
    figure
    subplot(2,2,1)
    imagesc(I_phase); axis image;
    subplot(2,2,2)
    imagesc(abs(f0).^0.2); axis image;
    
    subplot(2,2,3)
    imagesc(I_phase_new); axis image;
    subplot(2,2,4)
    imagesc(abs(f).^0.2); axis image; 
    colormap jet
end



end

