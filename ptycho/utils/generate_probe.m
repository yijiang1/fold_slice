function [probe] = generate_probe(N, lambda, dx, Ls, setup)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Parameters:    N      -> number of pixels
%               lambda -> the wave length 
%               dx     -> pixel size (in sample plane)
%               Ls     -> distance from focal plane to sample

    %Fresnel Zone Plate
    switch setup
        case "velo"
            Rn=90e-6; 
            dRn=50e-9;
        case "barry"
            Rn=80e-6; 
            dRn=70e-9;  
        case "barry2"
            Rn=70e-6;
            dRn=160e-9;
        otherwise
            Rn=90e-6;
            dRn=50e-9;
    end
    
    fl=2*Rn*dRn/lambda;%focal length corresponding to central wavelength
    D_FZP=180e-6;%dimeter of the FZP
    D_H=60e-6;%central beamstop

    %pixel size on FZP plane
    dx_fzp=lambda*fl/N/dx;
    %Coordinate on FZP plane
    lx_fzp=linspace(-dx_fzp*N/2,dx_fzp*N/2,N);
    [x_fzp,y_fzp]=meshgrid(lx_fzp);
    %Transmission function of the FZP
    T=exp(-1j*2*pi/lambda*(x_fzp.^2+y_fzp.^2)/2/fl); 
    C=double(sqrt(x_fzp.^2+y_fzp.^2)<=(D_FZP/2));% Cercular function of FZP
    H=double(sqrt(x_fzp.^2+y_fzp.^2)>=(D_H/2));%cental block

    %probe on sample plane
    probe=fresnel_propagation(C.*T.*H,dx_fzp,(fl+Ls),lambda);
    %figure(1);imagesc(abs(probe));axis image
    %figure(2);imagesc(angle(probe));axis image
end
