function [OUT]=fresnel_propagation(IN,dxy,z,lambda)
%FRESNEL_PROP Summary of this function goes here
%   Detailed explanation goes here
%Parameters:    dx,dy  -> the pixel pitch of the object
%               z      -> the distance of the propagation 
%               lambda -> the wave length 
%               X,Y    -> meshgrid of coordinate 
%               IN     -> input object

[M,N]=size(IN);
k=2*pi/lambda;

%coodinate grid
lx=linspace(-dxy*M/2,dxy*M/2,M);
[x,y]=meshgrid(lx);

% coordinate on the output plane
fc=1/dxy;
fu=lambda*z*fc;
lu=ifftshift(ifftshift(linspace(-fu/2,fu/2,M),1),2);

[u,v]=meshgrid(lu);


if z>0      
    pf=exp(1j*k*z)*exp(1j*k*(u.^2+v.^2)/2/z);
    kern=IN.*exp(1j*k*(x.^2+y.^2)/2/z);
    
    kerntemp=fftshift(kern);
    cgh=fft2(kerntemp);
    OUT=fftshift(fftshift(cgh.*pf,1),2);
else
    z=abs(z);
    pf = exp(1j*k*z)*exp(1j*k*(x.^2+y.^2)/2/z);
    cgh = ifft2(ifftshift(IN)./exp(1j*k*(u.^2+v.^2)/2/z));
    OUT = fftshift(cgh)./pf;      


end

