%view_probe_propagation_electron.m
% Visualize STEM probe propagation in free space
addpath(fullfile(pwd,'utils'))

%% specify parameters
N = 128; %probe/cbed size
voltage = 300; %beam voltage in keV
alpha_max = 30; %semi-convergence angle in mrad
dx = 0.1; %pixel size in angstrom
defocus = 100; %probe defocus in angstrom. NOTE: positive value gives underfocused probe
C3 = 1e5; %third-order spherical aberration in angstrom

z_distance = 10; %distance (in angstrom) between layers 
N_layers = 15; %number of layers

%%
close all
% generate inital probe
par_probe = {};
par_probe.df = defocus;
par_probe.C3 = C3;
par_probe.voltage = voltage;
par_probe.alpha_max = alpha_max;
par_probe.plotting = true;

[probe_init, ~] = make_tem_probe(dx, N, par_probe);

% generate free-space propagator
lambda = 12.398/sqrt((2*511.0+par_probe.voltage).*par_probe.voltage); %angstrom

[~,H,~,~] = near_field_evolution(ones(N), z_distance, lambda,  dx.*[N,N], true );
H = ifftshift(H);

% propagate probes
probes = zeros(N, N, N_layers);
probes(:,:,1) = probe_init;
for k=2:N_layers
    probes(:,:,k) = ifft2(bsxfun(@times, H, fft2(probes(:,:,k-1))));
end

% visualize probe magnitudes
imagesc3D(abs(probes)); axis image
