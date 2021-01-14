function [probe,mask,A] = generateProbeFunction(dx,N, px,py,df ,cs, N_z, Voltage, alpha_max, showFigure )
%Generate probe function
%   dx: pixel size in angstrom
%   N: number of pixels in each side
%   df: defocus in angstrom
%   Nz:
%   Voltage: beam voltage in keV
%   alpha_max: convergence angle in mrad


%Parameter
C3 = cs; %angstrom
C5 = 0;
C7 = 0;
lambda = 12.398/sqrt((2*511.0+Voltage).*Voltage); %angstrom

amax = alpha_max*1e-3; %in rad
amin = 0;

klimitmax = amax/lambda;
klimitmin = amin/lambda;
dk = 1/(dx*N);

if N_z==1
    dff = df;
else
    dff = linspace(-floor(N_z/2),ceil(N_z/2)-1,N_z)*df; %angstrom
end

kx = linspace(-floor(N/2),ceil(N/2)-1,N);
[kX,kY] = meshgrid(kx,kx);

kX = kX.*dk;
kY = kY.*dk;
kR = sqrt(kX.^2+kY.^2);

chi = zeros(N,N,N_z);
probe = zeros(N,N,N_z);

mask = single(kR<=klimitmax).*single(kR>=klimitmin);
%figure
%imshow(mask)
for i= 1:N_z
    chi(:,:,i) = -pi*lambda*kR.^2*dff(i) + pi/2*C3*lambda^3*kR.^4+pi/3*C5*lambda^5*kR.^6+pi/4*C7*lambda^7*kR.^8;
    
    phase = exp(-1i.*chi(:,:,i)).*exp(-2*pi*1i*px.*kX).*exp(-2*pi*1i*py.*kY);
    A = mask.*phase;
    probe(:,:,i) = mask.*phase;
    
    probe(:,:,i) = fftshift(ifft2(ifftshift(probe(:,:,i))));
    
    probe(:,:,i) = probe(:,:,i)/sum(sum(abs(probe(:,:,i))));
    if showFigure
        figure
        imagesc(abs(probe))
        axis image
    end
end

if showFigure
    figure
    imagesc(abs(mask))
    axis image

end
%probe = probe/(L/N)^2; %normalize the probe. Because in reciprocal space
%the zero frequency amplitute is 1, so, the summation is also one.
end

