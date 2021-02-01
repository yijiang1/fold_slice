function [ output] = shift( input,dx_x,dx_y,px,py )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


Ny = size(input,1);
Nx = size(input,2);

dk_y = 1/(dx_y*Ny);
dk_x = 1/(dx_x*Nx);

ky = linspace(-floor(Ny/2),ceil(Ny/2)-1,Ny);
kx = linspace(-floor(Nx/2),ceil(Nx/2)-1,Nx);

[kX,kY] = meshgrid(kx,ky);
kX = kX.*dk_x;
kY = kY.*dk_y;

f = fftshift(fft2(ifftshift(input)));
f = f.*exp(-2*pi*1i*px*kX).*exp(-2*pi*1i*py*kY);
output = fftshift(ifft2(ifftshift(f)));

end