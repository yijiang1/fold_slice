function [probe, mask] = make_tem_probe(dx, N, param)
%MAKE_TEM_PROBE Generate probe functions produced by object lens in 
% transmission electron microscope.
% Written by Yi Jiang dased on Eq.(2.10) in Advanced Computing in Electron 
% Microscopy (2nd edition) by Dr.Kirkland
% 
% Outputs:
%   probe: complex probe functions
%   mask: 2D Fourier mask function determined by alpha_max
%
% Inputs: 
%   dx: pixel size in angstrom
%   N: number of pixels in each side of the 2D probe
%   param: probe parameters and other settings

% parse inputs
parse_param = inputParser;
parse_param.KeepUnmatched = true;
parse_param.addParameter('df', 0, @isnumeric) %first-order aberration (defocus) in angstrom
parse_param.addParameter('C3', 0, @isnumeric) %third-order spherical aberration in angstrom
parse_param.addParameter('C5', 0, @isnumeric) %fifth-order spherical aberration in angstrom
parse_param.addParameter('C7', 0, @isnumeric) %seventh-order spherical aberration in angstrom
parse_param.addParameter('f_a2', 0, @isnumeric) %twofold astigmatism in angstrom
parse_param.addParameter('theta_a2', 0, @isnumeric) %azimuthal orientation in radian
parse_param.addParameter('f_a3', 0, @isnumeric) %threefold astigmatism in angstrom
parse_param.addParameter('theta_a3', 0, @isnumeric) %azimuthal orientation in radian
parse_param.addParameter('f_c3', 0, @isnumeric) %coma in angstrom
parse_param.addParameter('theta_c3', 0, @isnumeric) %azimuthal orientation in radian

parse_param.addParameter('shifts', [0,0], @isnumeric) %shift probe center in angstrom
parse_param.addParameter('N_z', 1, @isnumeric) %number of focal planes
parse_param.addParameter('voltage', 300, @isnumeric) %beam voltage
parse_param.addParameter('alpha_max', 30, @isnumeric) %semi-convergence angle
parse_param.addParameter('plotting', 0, @islogical) %plot probes

parse_param.parse(param)
p = parse_param.Results;

lambda = 12.398/sqrt((2*511.0+p.voltage).*p.voltage); %angstrom

amax = p.alpha_max*1e-3; %in rad
amin = 0;

klimitmax = amax/lambda;
klimitmin = amin/lambda;
dk = 1/(dx*N);

if p.N_z==1
    dff = p.df;
else
    dff = linspace(-floor(p.N_z/2),ceil(p.N_z/2)-1,p.N_z)*p.df; %angstrom
end

kx = linspace(-floor(N/2),ceil(N/2)-1,N);
[kX,kY] = meshgrid(kx,kx);

kX = kX.*dk;
kY = kY.*dk;
kR = sqrt(kX.^2+kY.^2);
theta = atan2(kY,kX);
chi = zeros(N,N,p.N_z);
probe = zeros(N,N,p.N_z);

mask = single(kR<=klimitmax).*single(kR>=klimitmin);
for i= 1:p.N_z
    chi(:,:,i) = -pi*lambda*kR.^2*dff(i);
    if p.C3~=0; chi(:,:,i) = chi(:,:,i) + pi/2*p.C3*lambda^3*kR.^4; end
    if p.C5~=0; chi(:,:,i) = chi(:,:,i) + pi/3*p.C5*lambda^5*kR.^6; end
    if p.C7~=0; chi(:,:,i) = chi(:,:,i) + pi/4*p.C7*lambda^7*kR.^8; end
    if p.f_a2~=0; chi(:,:,i) = chi(:,:,i) + pi*p.f_a2*lambda*kR.^2*sin(2*(theta-p.theta_a2)); end
    if p.f_a3~=0; chi(:,:,i) = chi(:,:,i) + 2*pi/3*p.f_a3*lambda^2*kR.^3*sin(3*(theta-p.theta_a3)); end
    if p.f_c3~=0; chi(:,:,i) = chi(:,:,i) + 2*pi/3*p.f_c3*lambda^2*kR.^3*sin(theta-p.theta_c3); end

    phase = exp(-1i.*chi(:,:,i)).*exp(-2*pi*1i*p.shifts(1).*kX).*exp(-2*pi*1i*p.shifts(2).*kY);
    probe(:,:,i) = mask.*phase;
    probe(:,:,i) = fftshift(ifft2(ifftshift(probe(:,:,i))));
    probe(:,:,i) = probe(:,:,i)/sum(sum(abs(probe(:,:,i))));
end

if p.plotting
    figure
    subplot(1,2,1)
    imagesc(abs(mask))
    axis image
    title('Fourier mask magnitude')

    subplot(1,2,2)
    colormap parula
    iaxis{1} = ([1 N]-floor(N/2)+1)*dx/10;
    iaxis{2} = ([1 N]-floor(N/2)+1)*dx/10;
    imagesc3D(iaxis{1},iaxis{2},abs(probe)) 
    xlabel('nm')
    ylabel('nm')
    axis image
      
    title('Probe magnitude')
end

%probe = probe/(L/N)^2; %normalize the probe. Because in reciprocal space
%the zero frequency amplitute is 1, so, the summation is also one.
end

