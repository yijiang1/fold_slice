% multislice_dof_summary.m
% summarize depth of focus reported in x-ray multislice ptycho literature

addpath(fullfile(pwd,'utils'))
dof = [];
pixel_size = [];
thickness = [];
dz = [];
Nlayer = [];
label = {};
alpha = 4;     % additional scaling coefficient

%% 1.1 Multi‐slice ptychography enables high‐resolution measurements in extended chemical reactors
% https://doi.org/10.1038/s41598-020-80926-6
energy = 9.1; %kev
det_pixel = 75e-6;  %detector pixel size (m) 
det_N = 512; % number of pixels in the detector
distance = 4.16; % # sample to detector distance (m)

[dof(1),pixel_size(1)] = DOF_calculator(energy, det_pixel, det_N, distance, alpha);
thickness(1) = 100e-6;
dz(1) = 100e-6; % layer distance used in multi-slice reconstruction
Nlayer(1) = 2; % number of layers used in multi-slice reconstruction
label{1} = 'PETRA III-polyimide foil';
%% 1.2 Multi‐slice ptychography enables high‐resolution measurements in extended chemical reactors
% https://doi.org/10.1038/s41598-020-80926-6
energy = 15.25; %kev
det_pixel = 75e-6;  %detector pixel size (m) 
det_N = 512; % number of pixels in the detector
distance = 3.435; % # sample to detector distance (m)

[dof(end+1),pixel_size(end+1)] = DOF_calculator(energy, det_pixel, det_N, distance, alpha);
thickness(end+1) = 650e-6;
dz(end+1) = 650e-6; % layer distance used in multi-slice reconstruction
Nlayer(end+1) = 2; % number of layers used in multi-slice reconstruction
label{end+1} = 'PETRA III-chemical reactor';
%% 2. Multi-slice ptychography with large numerical aperture multilayer Laue lenses
% https://doi.org/10.1364/OPTICA.5.000601
energy = 12; %kev
det_pixel = 55e-6;  %detector pixel size (m) 
det_N = 128; % number of pixels in the detector
distance = 0.5; % # sample to detector distance (m)

[dof(end+1),pixel_size(end+1)] = DOF_calculator(energy, det_pixel, det_N, distance, alpha);
thickness(end+1) = 10e-6;
dz(end+1) = 10e-6; % layer distance used in multi-slice reconstruction
Nlayer(end+1) = 2; % number of layers used in multi-slice reconstruction
label{end+1} = 'BNL';
%% 3.1 X-ray ptychography with extended depth of field - real data
% https://doi.org/10.1364/OE.24.029089
energy = 6.2; %kev
det_pixel = 172e-6;  %detector pixel size (m) 
det_N = 192; % number of pixels in the detector
distance = 7.2; % # sample to detector distance (m)

[dof(end+1),pixel_size(end+1)] = DOF_calculator(energy, det_pixel, det_N, distance, alpha);
thickness(end+1) = 200e-6;
dz(end+1) = 200e-6; % layer distance used in multi-slice reconstruction
Nlayer(end+1) = 2; % number of layers used in multi-slice reconstruction
label{end+1} = 'PSI';
%% 3.2 X-ray ptychography with extended depth of field - simulation
% https://doi.org/10.1364/OE.24.029089
energy = 6.2; %kev
det_pixel = 172e-6;  %detector pixel size (m) 
det_N = 512; % number of pixels in the detector
distance = 7.2; % # sample to detector distance (m)

[dof(end+1),pixel_size(end+1)] = DOF_calculator(energy, det_pixel, det_N, distance, alpha);
thickness(end+1) = [40e-6];
dz(end+1) = [20e-6]; % layer distance used in multi-slice reconstruction
Nlayer(end+1) = 3; % number of layers used in multi-slice reconstruction
label{end+1} = 'PSI-sim';
%% 4. High-Resolution Multislice X-Ray Ptychography of Extended Thick Objects
% https://doi.org/10.1103/PhysRevLett.112.053903
energy = 7; %kev
dx = 20e-9;
det_pixel = 75e-6;  %detector pixel size (m). Not given in the paper
det_N = 606; % number of pixels in the detector
distance = dx*det_pixel*det_N/(1.23984193e-9/energy); % # sample to detector distance (m)

[dof(end+1),pixel_size(end+1)] = DOF_calculator(energy, det_pixel, det_N, distance, alpha);
thickness(end+1) = 105e-6;
dz(end+1) = 1e-9; % layer distance used in multi-slice reconstruction
Nlayer(end+1) = 1; % number of layers used in multi-slice reconstruction
label{end+1} = 'SPring 8';
%% 5. Resolving 500 nm axial separation by multi-slice X-ray ptychography
% https://doi.org/10.1107/S2053273318017229
energy = 12; %kev
det_pixel = 55e-6;  %detector pixel size (m)
det_N = 300; % number of pixels in the detector
distance = 0.35; % # sample to detector distance (m)

[dof(end+1),pixel_size(end+1)] = DOF_calculator(energy, det_pixel, det_N, distance, alpha);
thickness(end+1) = 500e-9;
dz(end+1:end+1) = 500e-9; % layer distance used in multi-slice reconstruction
Nlayer(end+1) = 2; % number of layers used in multi-slice reconstruction
label{end+1} = 'BNL-XRF';
%% 6. 3D x-ray imaging of continuous objects beyond the depth of focus limit - simulation, tomography
% https://doi.org/10.1364/OE.24.029089
energy = 5; %kev
dx = 1e-9;
det_pixel = 75e-6;  %detector pixel size (m). Not given in the paper 
det_N = 72; % number of pixels in the detector
distance = dx*det_pixel*det_N/(1.23984193e-9/energy); % # sample to detector distance (m)

[dof(end+1),pixel_size(end+1)] = DOF_calculator(energy, det_pixel, det_N, distance, alpha);
thickness(end+1) = 160e-9;
dz(end+1:end+1) = 1e-9; % layer distance used in multi-slice reconstruction
Nlayer(end+1) = 1; % number of layers used in multi-slice reconstruction
label{end+1} = 'sim-tomo';

%% 7. Adorym: A multi-platform generic x-ray image reconstruction framework based on automatic differentiation
%%%%%%% I'm not sure if they really used multi-slice...
%{
% https://arxiv.org/abs/2012.12686
energy = 5.5; %kev
det_pixel = 172e-6;  %detector pixel size (m). Not given in the paper 
det_N = 64; % number of pixels in the detector
distance = 2; % # sample to detector distance (m)

[dof(end+1),pixel_size(end+1)] = DOF_calculator(energy, det_pixel, det_N, distance, alpha);
thickness(end+1) = 8e-6;
dz(end+1:end+1) = 1e-9; % layer distance used in multi-slice reconstruction
Nlayer(end+1) = 1; % number of layers used in multi-slice reconstruction
label{end+1} = 'APS-bnp-algae???';
%}
%% LCO
% 
energy = 9.3; %kev
det_pixel = 75e-6;  %detector pixel size (m) 
det_N = 64; % number of pixels in the detector
distance = 1.92; % # sample to detector distance (m)
alpha = 4;     % additional scaling coefficient

[dof(end+1),pixel_size(end+1)] = DOF_calculator(energy, det_pixel, det_N, distance, alpha);
thickness(end+1) = 25e-6;
dz(end+1) = 8e-6; % layer distance used in multi-slice reconstruction
Nlayer(end+1) = 5; % number of layers used in multi-slice reconstruction
label{end+1} = 'APS-velo-LCO';

%% IC Pillar
% 
energy = 8.8; %kev
det_pixel = 75e-6;  %detector pixel size (m) 
det_N = 64; % number of pixels in the detector
distance = 1.92; % # sample to detector distance (m)
alpha = 4;     % additional scaling coefficient

[dof(end+1),pixel_size(end+1)] = DOF_calculator(energy, det_pixel, det_N, distance, alpha);
thickness(end+1) = 20e-6;
dz(end+1) = 8e-6; % layer distance used in multi-slice reconstruction
Nlayer(end+1) = 5; % number of layers used in multi-slice reconstruction
label{end+1} = 'APS-velo-IC-pillar';

%% CIGS
% 
energy = 9.3; %kev
det_pixel = 75e-6;  %detector pixel size (m) 
det_N = 64; % number of pixels in the detector
distance = 1.92; % # sample to detector distance (m)
alpha = 4;     % additional scaling coefficient

[dof(end+1),pixel_size(end+1)] = DOF_calculator(energy, det_pixel, det_N, distance, alpha);
thickness(end+1) = 30e-6;
dz(end+1) = 8e-6; % layer distance used in multi-slice reconstruction
Nlayer(end+1) = 5; % number of layers used in multi-slice reconstruction
label{end+1} = 'APS-velo-CIGS';

%%
close all
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
for i=1:length(thickness)
    plot(i,thickness(i)./dof(i),'.','MarkerSize',15, 'DisplayName',label{i})
end
%hold on
line([1 length(dof)],[1 1],'LineWidth',1,'LineStyle','--','Color','r', 'DisplayName', 'DOF')
legend
set(axes1,'YMinorTick','on','YScale','log');
