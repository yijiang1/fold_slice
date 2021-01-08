 

%load('/home/beams0/YJIANG/ptychography/electron/Nb3Cl8/0205_83/data/data_Nb3Cl8_0205_83_roi3.mat');
%load('/home/beams0/YJIANG/ptychography/electron/Nb3Cl8/19/data/data_Nb3Cl8_19_roi8.mat');
load('//home/beams0/YJIANG/ptychography/electron/PrScO3/18/data/data_PSO_18_roi3.mat');

%load('//home/beams0/YJIANG/ptychography/electron/mos2/22/data_mos2_22_roi6_pos.mat');
roi = '3';
Ndp = 256;
rot_angle = 30;

bg_level = 0;
transpose = true;

%%
dp(dp<bg_level) = 0;
Ny = size(dp,3);
Nx = size(dp,4);
dp = reshape(dp,size(dp,1), size(dp,2),Ny*Nx);

if transpose
    dp = permute(dp,[2,1,3]);
end
if Ndp~=size(dp,1)
    dp = crop_pad_3D( dp, [Ndp,Ndp,size(dp,3)]);
end

saveName = strcat('data_roi',roi,'_dp.hdf5');
h5create(saveName, '/dp', size(dp),'ChunkSize',[Ndp Ndp min([Ny,Nx,100])],'Deflate',4)
h5write(saveName, '/dp', dp)
%{
py = linspace(1,Ny,Ny)*scanStepSize_y;
%py = py - mean(py);
px = linspace(1,Nx,Nx)*scanStepSize_x;
%px = px - mean(px);

[ppX0,ppY0] = meshgrid(px,py);
ppY_rot = ppX0*-sind(rot_angle) + ppY0*cosd(rot_angle);
ppX_rot = ppX0*cosd(rot_angle) + ppY0*sind(rot_angle);
ppX_rot = ppX_rot(:);
ppY_rot = ppY_rot(:);
%ppX_rot = ppX_rot - mean(ppX_rot);
%ppY_rot = ppY_rot - mean(ppY_rot);

saveName = strcat('data_roi',roi,'_para.hdf5');

hdf5write(saveName, '/ppX', ppX(:))
hdf5write(saveName, '/ppY', ppY(:),'WriteMode','append')    
%}