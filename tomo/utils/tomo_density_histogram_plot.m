% tomo_quantitative.m
import plotting.franzmap

matlab_tomo_path='/mnt/das-gpfs/work/p16167/matlab_new/tomo/';
cd(matlab_tomo_path)
addpath([matlab_tomo_path 'utils'])
return


%% Constants:

scrsz = get(0,'ScreenSize');
tomo_folder='tomo_S03041_to_S04042_500x500_run_1_c';
tomo_path_read= ['/sls/X12SA/Data20/e16167/analysis_tomo/' tomo_folder '/']; 
tomo_file_name='tomogram_delta_S03041_S04042_Hann_freqscl_1.00.mat';

%% Read tomographic reconstruction:

load([tomo_path_read tomo_file_name]);
tomo_path_write= sprintf('/mnt/das-gpfs/work/p16167/analysis_tomo_offline/%s/quantitative_%s_%4.2f/',tomo_folder,filter_type,freq_scale); 

%% Make histogram of whole sample

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Choose parameters
sam=1000;      % Number of bins in histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove data ouside computed tomogram
N = size(tomogram_delta,1);
xt = [-N/2:N/2-1];
[Xt Yt] = meshgrid(xt,xt);
circulo = 1-radtap(Xt,Yt,10,N/2-3);
cylinder=repmat(circulo,[1 1 size(tomogram_delta,3)]);
data=tomogram_delta.*cylinder;

% Calculate whole histogram
M=size(data,1)*size(data,2)*size(data,3);
data_long=reshape(data,M,1);
cylinder_long=reshape(cylinder,M,1);
data_nozeros=data_long(cylinder_long == 1);
[hst,bins]=hist(data_nozeros,sam);

clear circulo
clear cylinder
clear tomogram_delta

%% Plot histogram

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Choose parameters
quant='eden';  % Choose quantity to plot: 'delta' for delta or 'eden' for electron density
yaxis='log';   % Y axis can be linear ('lin') or logaritmic ('log')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);

if isstr(quant)&&strcmpi('eden',quant)
    bins_plot=bins*factor_edensity;
    xaxis_label='electron density (A^{-3})';
elseif isstr(quant)&&strcmpi('delta',quant)
    bins_plot=bins;
    xaxis_label='delta';
else
    error('Supported strings for quant are delta or eden')
end

if isstr(yaxis)&&strcmpi('lin',yaxis)
    plot(bins_plot,hst); xlabel(xaxis_label); ylabel('number of voxels');
elseif  isstr(yaxis)&&strcmpi('log',yaxis)
    semilogy(bins_plot,hst); xlabel(xaxis_label); ylabel('number of voxels');
else
    error('Supported strings for yaxis are lin or log')
end 

%% Save histogram data

savedata=0;  % Equal to 1 for saving data, or 0 for not saving

if savedata
    fid=fopen([tomo_path_write sprintf('histogram_%s.txt',tomo_folder)],'w');
    fprintf(fid, '# delta \t electron density (Angtrom-3) \t  number of voxels\n');
    for hh=1:length(bins)
        fprintf(fid, '%e \t %e \t %e\n', bins(hh),factor_edensity*bins(hh),hst(hh));
    end
    fclose(fid)
    
    save(sprintf('%shistogram_%s.mat',tomo_path_write,tomo_folder),'bins','factor_edensity','hst','tomo_path_read');
   
    print('-f1','-depsc2', [ tomo_path_write sprintf('histogram_%s.eps',tomo_folder)]);
    print('-f1','-dpng', [ tomo_path_write sprintf('histogram_%s.png',tomo_folder)]);

end
    

%% Plot slices to navigate in 3D data (with color lines)

%Choose parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
analysis_case='cell1_nucleolus';   % please chose a different name for different slected volumes to save data in separate folders
quant='eden';             % choose quantity to plot: 'delta' for delta or 'eden' for electron density
scl=[0.25 0.45];               % color scale can be 'auto' for automatic or e.g. [0.25 0.45]
valz=80;                  % z coordinate to select slice in xy plane
valx=797;                 % x coordinate to select slice in yz plane
valy=795;                 % y coordinate to select slice in xz plane
sidex=20;                 % box size in x for volume of interest
sidey=20;                 % box size in y for volume of interest
sidez=20;                 % box size in z for volume of interest
colorx='r';
colory='b';
colorz='g';
color_map='jet';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xs=valx-round(sidex/2);
xf=valx+round(sidex/2);
ys=valy-round(sidey/2);
yf=valy+round(sidey/2);
zs=valz-round(sidez/2);
zf=valz+round(sidez/2);

if isstr(quant)&&strcmpi('eden',quant)
    data_corr=data*factor_edensity;
    xaxis_label='electron density (A^{-3})';
elseif isstr(quant)&&strcmpi('delta',quant)
    data_corr=data;
    xaxis_label='delta';
else
    error('Supported strings for quant are delta or eden')
end

if isstr(scl)&&strcmpi('auto',scl)
    scale=[min(data_corr(:)) max(data_corr(:))];
else
    scale=scl;
end

figure(2);
%figure('Position',[1,400,800,800]); 
subplot(2,2,3);
imagesc(data_corr(:,:,valz), scale); axis xy equal tight;
xlabel('x'); ylabel('y')
title(sprintf('z = %d',valz)); colormap bone(256); hold on;
plot([valx,valx],[1,size(data_corr,1)],colorx);
plot([1,size(data_corr,2)],[valy,valy],colory);
plot([1,size(data_corr,2)],[1,1],colorz,'Linewidth',3);
plot([1,size(data_corr,2)],[size(data_corr,1),size(data_corr,1)],colorz,'Linewidth',3);
plot([1,1],[1,size(data_corr,1)],colorz,'Linewidth',3);
plot([size(data_corr,2),size(data_corr,2)],[1,size(data_corr,1)],colorz,'Linewidth',3);
plot([xs,xf],[ys,ys],colorz);
plot([xs,xf],[yf,yf],colorz);
plot([xs,xs],[ys,yf],colorz);
plot([xf,xf],[ys,yf],colorz);
hold off;

subplot(2,2,4);
imageyz=(squeeze(data_corr(:,valx,:)));
imagesc(imageyz, scale); axis xy equal tight; colorbar;
xlabel('z'); ylabel('y');
title(sprintf('x = %d',valx)); colormap bone(256); hold on;
plot([valz,valz],[1,size(data_corr,1)],colorz);
plot([1,size(data_corr,3)],[valy,valy],colory);
plot([1,size(data_corr,3)],[1,1],colorx,'Linewidth',3);
plot([1,size(data_corr,3)],[size(data_corr,1),size(data_corr,1)],colorx,'Linewidth',3);
plot([1,1],[1,size(data_corr,1)],colorx,'Linewidth',3);
plot([size(data_corr,3),size(data_corr,3)],[1,size(data_corr,1)],colorx,'Linewidth',3);
plot([zs,zf],[ys,ys],colorx);
plot([zs,zf],[yf,yf],colorx);
plot([zs,zs],[ys,yf],colorx);
plot([zf,zf],[ys,yf],colorx);
hold off;

subplot(2,2,1);
imagexz=(squeeze(data_corr(valy,:,:)))';
imagesc(imagexz, scale); axis xy equal tight;
xlabel('x'); ylabel('z');
title(sprintf('y = %d',valy)); colormap bone(256); hold on;
plot([valx,valx],[1,size(data_corr,3)],colorx);
plot([1,size(data_corr,2)],[valz,valz],colorz);
plot([1,size(data_corr,2)],[1,1],colory,'Linewidth',3);
plot([1,size(data_corr,2)],[size(data_corr,3),size(data_corr,3)],colory,'Linewidth',3);
plot([1,1],[1,size(data_corr,3)],colory,'Linewidth',3);
plot([size(data_corr,2),size(data_corr,2)],[1,size(data_corr,3)],colory,'Linewidth',3);
plot([xs,xf],[zs,zs],colory);
plot([xs,xf],[zf,zf],colory);
plot([xs,xs],[zs,zf],colory);
plot([xf,xf],[zs,zf],colory);
hold off;
set(gcf,'Outerposition',[1 5 800 800])

%% Histogram of selected voi

% Choose parameters %%%%%%%%%%%%%%%%
sampling_sel=50;              % Number of bins in histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_sel=data_corr(ys:yf,xs:xf,zs:zf);

figure(3);
%figure('Position',[1,400,800,800]); 
subplot(2,2,3);
imagesc(data_sel(:,:,round((zf-zs)/2)), scale); axis xy equal tight;
xlabel('x'); ylabel('y')
title(sprintf('z = %d',valz)); colormap bone(256); hold on;
hold off;

subplot(2,2,4);
imagesc(squeeze(data_sel(:,round((xf-xs)/2),:)), scale); axis xy equal tight;
xlabel('z'); ylabel('y');
title(sprintf('x = %d',valx)); colormap bone(256); hold on;
hold off;

subplot(2,2,1);
imagesc(squeeze(data_sel(round((yf-ys)/2),:,:))', scale); axis xy equal tight;
xlabel('x'); ylabel('z');
title(sprintf('y = %d',valy)); colormap bone(256); hold on;
hold off;

M_sel=size(data_sel,1)*size(data_sel,2)*size(data_sel,3);
data_sel_long=reshape(data_sel,M_sel,1);
[hst_sel,bins_sel]=hist(data_sel_long,sampling_sel);

figure(4)
plot(bins_sel, hst_sel)
xlabel(xaxis_label)
ylabel('number of voxels')
title('histogram of VOI')

%% Make individual plots without lines

x=((1:size(data_corr,2))-round(size(data_corr,2))/2)*pixsize*1e6; % [microns]
y=((1:size(data_corr,1))-round(size(data_corr,1))/2)*pixsize*1e6; % [microns]
z=((1:size(data_corr,3))-round(size(data_corr,3))/2)*pixsize*1e6; % [microns]

figure(5)
imagexz=(squeeze(data_corr(valy,:,:)))';
imagesc(x,z,imagexz, scale); axis xy equal tight;
xlabel('x (microns)'); ylabel('z (microns)');
title(sprintf('electron density [e/A^3]; y = %d',valy)); colormap bone(256); 
colorbar;

figure(6)
imageyz=(squeeze(data_corr(:,valx,:)))';
imagesc(y,z,imageyz, scale); axis xy equal tight; colorbar;
xlabel('y (microns)'); ylabel('z (microns)');
title(sprintf('electron density [e/A^3]; x = %d',valx)); colormap bone(256); 
colorbar;

figure(7)
imagesc(x,y,data_corr(:,:,valz), scale); axis xy equal tight;
title(sprintf('electron density [e/A^3]; z = %d',valz)); colormap bone(256); 
xlabel('x (microns)'); ylabel('y (microns)');
colorbar

%% Fit histogram peak to Gaussian curve

fit_type='gauss2'; % try 'gauss1' for one peak and 'gauss2' for a double peak fit

f = fit(bins_sel.',hst_sel.',fit_type)  
figure(8)
plot(f,bins_sel,hst_sel)

value=f.b1;
sigma=f.c1/sqrt(2);
FWHM=2.35482*sigma;

if isstr(quant)&&strcmpi('eden',quant)
    display(sprintf('electron density: %f4.2 +/- %f4.2',value,sigma))
else isstr(quant)&&strcmpi('delta',quant)
    display(sprintf('delta: %e +/- %e',value*factor_edensity,sigma*factor_edensity))
end

if isstr(fit_type)&&strcmpi('gauss2',fit_type)
    value2=f.b2;
    sigma2=f.c2/sqrt(2);
    FWHM2=2.35482*sigma2;
    if isstr(quant)&&strcmpi('eden',quant)
        display(sprintf('electron density: %f4.2 +/- %f4.2',value2,sigma2))
    else isstr(quant)&&strcmpi('eden',quant)
        display(sprintf('delta: %e +/- %e',value2*factor_edensity,sigma2*factor_edensity))
    end
    
end
    
%% Estimate mass density
% Choose parameters %%%%%%%%%%%%%%%%%%%
AZ_ratio=1.85; % Estimation of molar mass (g/mol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NA=6.022e23;  %[mol-1]

if isstr(quant)&&strcmpi('eden',quant)
   mass_density=value*AZ_ratio/NA*1e24;
   mass_density_sigma=sigma*AZ_ratio/NA*1e24
else isstr(quant)&&strcmpi('delta',quant)
    mass_density=value*factor_edensity*AZ_ratio/NA*1e24;
    mass_density_sigma=sigma*factor_edensity*AZ_ratio/NA*1e24
end

if isstr(fit_type)&&strcmpi('gauss2',fit_type)
  if isstr(quant)&&strcmpi('eden',quant)
   mass_density2=value2*AZ_ratio/NA*1e24;
   mass_density_sigma2=sigma2*AZ_ratio/NA*1e24
else isstr(quant)&&strcmpi('delta',quant)
    mass_density2=value2*factor_edensity*AZ_ratio/NA*1e24;
    mass_density_sigma2=sigma2*factor_edensity*AZ_ratio/NA*1e24
  end
end

display(sprintf('mass density: %f4.2 +/- %f4.2',mass_density,mass_density_sigma))

if isstr(fit_type)&&strcmpi('gauss2',fit_type)
   display(sprintf('mass density: %f4.2 +/- %f4.2',mass_density2,mass_density_sigma2))
end
%% Save analysis
savedata=1;
casefolder=[tomo_path_write analysis_case '/'];
savename=['histogram_VOI_' analysis_case];
if savedata == 1
    if ~exist('casefolder','dir'); mkdir(casefolder); end
    print('-f2','-depsc2', [ casefolder savename '_3D_orientation_all.eps']);
    print('-f2','-dtiff', [ casefolder savename '_3D_orientation_all.tif']);
    print('-f3','-depsc2', [ casefolder savename '_3D_orientation.eps']);
    print('-f3','-dtiff', [ casefolder savename '_3D_orientation.tif']);
    print('-f4','-depsc2', [ casefolder savename '_histogram.eps']);
    print('-f4','-dtiff', [ casefolder savename '_histogram.tif']);
    print('-f5','-depsc2', [ casefolder savename '_slice_y.eps']);
    print('-f5','-dtiff', [ casefolder savename '_slice_y.tif']);
    print('-f6','-depsc2', [ casefolder savename '_slice_x.eps']);
    print('-f6','-dtiff', [ casefolder savename '_slice_x.tif']);
    print('-f7','-depsc2', [ casefolder savename '_slice_z.eps']);
    print('-f7','-dtiff', [ casefolder savename '_slice_z.tif']);
    print('-f8','-depsc2', [ casefolder savename '_Gauss_fit.eps']);
    print('-f8','-dtiff', [ casefolder savename '_Gauss_fit.tif']);

    fid=fopen([casefolder savename 'histogram.txt'],'w');
    fprintf(fid, '# electron density (Angtrom-3) /  number of voxels\n');
    for hh=1:length(bins_sel)
        fprintf(fid, '%e %e\n', bins_sel(hh),hst_sel(hh));
    end
    fclose(fid)
    
    save([casefolder savename '.m'],'bins_sel','hst_sel','valx','valy',...
          'valz','sidex','sidey','sidez','analysis_case','tomo_path_read',...
          'tomo_path_write','output_folder','pixsize','sampling_sel','scale',...
          'quant','f','value','sigma','FWHM','AZ_ratio','mass_density','mass_density_sigma');
    if isstr(fit_type)&&strcmpi('gauss2',fit_type)
        save([casefolder savename '.m'],'bins_sel','hst_sel','valx','valy',...
          'valz','sidex','sidey','sidez','analysis_case','tomo_path_read',...
          'tomo_path_write','output_folder','pixsize','sampling_sel','scale',...
          'quant','f','value','sigma','FWHM','AZ_ratio','mass_density','mass_density_sigma',...
          'value2','sigma2','FWHM2','mass_density2','mass_density_sigma2');
    end
      
end



   
%% Delete large variables
% After this the code needs to be run from the very beginning to read the
% full tomogram
clear data0
clear data_corr


% %% Read amplitude data:
% 
% % This needs to be changed for each sample:
% filename_amp=[tomo_path 'tomogram_beta_S04693_S05999_Hann_freqscl_0.35.mat']; % tomorec
% ampdata = load(filename_amp) ;
% data_amp_sel=ampdata.tomogram_beta(ys:yf,xs:xf,zs:zf);
% 
% %% Plot full amplitude slices to navigate in 3D data (with color lines)
% 
% scale_amp=[-0.1e-6,1.3e-6];
% 
% figure(11);
% %figure('Position',[1,400,800,800]); 
% subplot(2,2,3);
% imagesc(ampdata.tomogram_beta(:,:,valz), scale_amp); axis xy equal tight;
% xlabel('x'); ylabel('y')
% title(sprintf('z = %d',valz)); colormap bone(256); hold on;
% plot([valx,valx],[1,size(ampdata.tomogram_beta,1)],colorx);
% plot([1,size(ampdata.tomogram_beta,2)],[valy,valy],colory);
% plot([1,size(ampdata.tomogram_beta,2)],[1,1],colorz,'Linewidth',3);
% plot([1,size(ampdata.tomogram_beta,2)],[size(ampdata.tomogram_beta,1),size(ampdata.tomogram_beta,1)],colorz,'Linewidth',3);
% plot([1,1],[1,size(ampdata.tomogram_beta,1)],colorz,'Linewidth',3);
% plot([size(ampdata.tomogram_beta,2),size(ampdata.tomogram_beta,2)],[1,size(ampdata.tomogram_beta,1)],colorz,'Linewidth',3);
% plot([xs,xf],[ys,ys],colorz);
% plot([xs,xf],[yf,yf],colorz);
% plot([xs,xs],[ys,yf],colorz);
% plot([xf,xf],[ys,yf],colorz);
% hold off;
% 
% subplot(2,2,4);
% imageyz=(squeeze(ampdata.tomogram_beta(:,valx,:)));
% imagesc(imageyz, scale_amp); axis xy equal tight; colorbar;
% xlabel('z'); ylabel('y');
% title(sprintf('x = %d',valx)); colormap bone(256); hold on;
% plot([valz,valz],[1,size(ampdata.tomogram_beta,1)],colorz);
% plot([1,size(ampdata.tomogram_beta,3)],[valy,valy],colory);
% plot([1,size(ampdata.tomogram_beta,3)],[1,1],colorx,'Linewidth',3);
% plot([1,size(ampdata.tomogram_beta,3)],[size(ampdata.tomogram_beta,1),size(ampdata.tomogram_beta,1)],colorx,'Linewidth',3);
% plot([1,1],[1,size(ampdata.tomogram_beta,1)],colorx,'Linewidth',3);
% plot([size(ampdata.tomogram_beta,3),size(ampdata.tomogram_beta,3)],[1,size(ampdata.tomogram_beta,1)],colorx,'Linewidth',3);
% plot([zs,zf],[ys,ys],colorx);
% plot([zs,zf],[yf,yf],colorx);
% plot([zs,zs],[ys,yf],colorx);
% plot([zf,zf],[ys,yf],colorx);
% hold off;
% 
% subplot(2,2,1);
% imagexz=(squeeze(ampdata.tomogram_beta(valy,:,:)))';
% imagesc(imagexz, scale_amp); axis xy equal tight;
% xlabel('x'); ylabel('z');
% title(sprintf('y = %d',valy)); colormap bone(256); hold on;
% plot([valx,valx],[1,size(ampdata.tomogram_beta,3)],colorx);
% plot([1,size(ampdata.tomogram_beta,2)],[valz,valz],colorz);
% plot([1,size(ampdata.tomogram_beta,2)],[1,1],colory,'Linewidth',3);
% plot([1,size(ampdata.tomogram_beta,2)],[size(ampdata.tomogram_beta,3),size(ampdata.tomogram_beta,3)],colory,'Linewidth',3);
% plot([1,1],[1,size(ampdata.tomogram_beta,3)],colory,'Linewidth',3);
% plot([size(ampdata.tomogram_beta,2),size(ampdata.tomogram_beta,2)],[1,size(ampdata.tomogram_beta,3)],colory,'Linewidth',3);
% plot([xs,xf],[zs,zs],colory);
% plot([xs,xf],[zf,zf],colory);
% plot([xs,xs],[zs,zf],colory);
% plot([xf,xf],[zs,zf],colory);
% hold off;
% set(gcf,'Outerposition',[1 300 800 800])
% %% Make individual plots without lines
% 
% figure(12)
% imagexz=(squeeze(ampdata.tomogram_beta(valy,:,:)))';
% imagesc(x,z,imagexz, scale_amp); axis xy equal tight;
% xlabel('x (microns)'); ylabel('z (microns)');
% title(sprintf('electron density [e/A^3]; y = %d',valy)); colormap bone(256); 
% colorbar;
% 
% figure(13)
% imageyz=(squeeze(ampdata.tomogram_beta(:,valx,:)))';
% imagesc(y,z,imageyz, scale_amp); axis xy equal tight; colorbar;
% xlabel('y (microns)'); ylabel('z (microns)');
% title(sprintf('electron density [e/A^3]; x = %d',valx)); colormap bone(256); 
% colorbar;
% 
% figure(14)
% imagesc(x,y,ampdata.tomogram_beta(:,:,valz), scale_amp); axis xy equal tight;
% title(sprintf('electron density [e/A^3]; z = %d',valz)); colormap bone(256); 
% xlabel('x (microns)'); ylabel('y (microns)');
% colorbar
% 
% %% Histogram of selected amplitude voi
% sampling_amp_sel=70;
% scale_amp=[-0.1e-6,1.3e-6];
% 
% figure(8);
% %figure('Position',[1,400,800,800]); 
% subplot(2,2,3);
% imagesc(data_amp_sel(:,:,round((zf-zs)/2)), scale_amp); axis xy equal tight;
% xlabel('x'); ylabel('y')
% title(sprintf('z = %d',valz)); colormap bone(256); hold on;
% hold off;
% 
% subplot(2,2,4);
% imagesc(squeeze(data_amp_sel(:,round((xf-xs)/2),:)), scale_amp)
% xlabel('z'); ylabel('y');
% title(sprintf('x = %d',valx)); colormap bone(256); hold on;
% hold off;
% 
% subplot(2,2,1);
% imagesc(squeeze(data_amp_sel(round((yf-ys)/2),:,:))', scale_amp)
% xlabel('x'); ylabel('z');
% title(sprintf('y = %d',valy)); colormap bone(256); hold on;
% hold off;
% 
% M_amp_sel=size(data_amp_sel,1)*size(data_amp_sel,2)*size(data_amp_sel,3);
% data_amp_sel_long=reshape(data_amp_sel,M_amp_sel,1);
% [hst_amp_sel,bins_amp_sel]=hist(data_amp_sel_long,sampling_amp_sel);
% 
% figure(9)
% plot(bins_amp_sel, hst_amp_sel)
% xlabel('beta')
% ylabel('number of voxels')
% title('histogram of VOI')
% %% Add path for Franzmap
% addpath('/mnt/das-gpfs/work/p15232/matlab/');
% %% Make bivariate histogram of voi
% 
% bins = 256; % number of bins of the histogram
% spacing = 'lin'; %'lin'; 'log'; % lin is better
% 
% delta_slice=data_sel./factor_edensity;
% abs_slice=data_amp_sel;
% 
% % find indices corresponding to the materials phase only (exclude air)
% % clear mask mask_ind
% mask=data_sel>1E-6;
% mask_ind=find(delta_slice>1E-6);
% 
% % Reshape the images into 1D vectors
% x=abs_slice(mask_ind);
% y=delta_slice(mask_ind);
% 
% clear xedges yedges
% switch lower(spacing)
%     case 'lin'
%         % linearly spaced edges of the histogram
%         xedges = linspace(min(x),max(x)+0.11e-6,bins);
%         yedges = linspace(min(y),max(y),bins);
%     case 'log'
%         xedges = linspace(min(x),max(x),bins);
%         yedges = logspace(log10(min(y)),log10(max(y)),bins);
% end
% 
% % Calculate the 1D histogram
% [xn, xbin] = histc(x,xedges);
% [yn, ybin] = histc(y,yedges);
% 
% %xbin, ybin zero for out of range values 
% % (see the help of histc) force this event to the 
% % first bins
% xbin(find(xbin == 0)) = inf;
% ybin(find(ybin == 0)) = inf;
% 
% xnbin = length(xedges);
% ynbin = length(yedges);
% 
% if xnbin >= ynbin
%     xy = ybin*(xnbin) + xbin;
%       indexshift =  xnbin; 
% else
%     xy = xbin*(ynbin) + ybin;
%       indexshift =  ynbin; 
% end
% 
% %[xyuni, m, n] = unique(xy);
% xyuni = unique(xy);
% xyuni(end) = []; 
% hstres = histc(xy,xyuni);
% clear xy;
% 
% histmat = zeros(ynbin,xnbin);
% histmat(xyuni-indexshift) = hstres;
% % %% Add path for Franzmap
% addpath('/afs/psi.ch/project/cxs/matlab/cSAXS_matlab_base_package/');
% %% display the bivariate histogram
% figure(10)
% sub1=subplot(3,3,[4,5,7,8]);
% imagesc(xedges.*1e7, yedges.*1e5, log10(histmat')), axis xy square tight
% xlim([0 14]);
% ylim([0.1 2.2]);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% For drawing the lines %%%%
% %lineh1= 1.0; % value in beta .*1e-5
% %lineh2= 1.0; % value in beta .*1e-5
% %linev1= 2.5; % value in delta .*1e-7
% %linev2= 2.5; % value in delta .*1e-7
% %%%%%%%   end of edit %%%%%%%%%%%%%%%%
% hold on
% %plot([-4 12],[lineh1 lineh1],'-b')
% %plot([-4 12],[lineh2 lineh2],'-b')
% %plot([linev1 linev1],[0.2 2],'-b')
% %plot([linev2 linev2],[0.2 2],'-b')
% hold off
% 
% thisfontsize=12;
% 
% colormap('franzmap')
% Contours =[1e0 1e1 1e2 1e3 1e4 1e5 1e6 1e7];
% hColorbar = colorbar('East','YTick',log10(Contours),'YTickLabel',Contours);
% hXLabel = xlabel('Absorption index, \beta [x 10^{-7}]');
% hYLabel = ylabel('Refractive index decrement, \delta [x 10^{-5}] ');
% set(gca,...
%     'FontName'  , 'Helvetica',...
%     'FontSize'  , thisfontsize         ,...
%     'Box'       , 'off'      ,...
%     'OuterPosition', [0 0 0.53 0.73] ,...
%     'TickDir'   , 'out'      ,...'YAxisLocation','right'
%     'XMinorTick', 'on'       ,...
%     'YMinorTick', 'on'       ,...
%     'XColor'    , [.0 .0 .0] ,...   
%     'YColor'    , [.0 .0 .0] ,...
%     'YTick'     ,  0:0.2:2.2   ,...
%     'XTick'     , -6:2:20    ,...
%     'LineWidth' , 1         );
% set([hXLabel,hYLabel],...
%     'FontName', 'Arial',...
%     'FontSize', thisfontsize-1    );
% set(hColorbar,...
%     'Box'    , 'on'         ,...
%     'TickDir', 'in'         ,...
%     'Direction','normal', ...
%     'YAxisLocation','left',...
%     'YColor' , [0.9 0.9 0.9] ,...
%     'XColor' , [0.9 0.9 0.9] , ...    
%     'Position',[0.47 0.11 0.03 0.3]);
% 
% subplot(3,3,[1,2])
% b=bar(xedges.*1e7,xn*.1e-5,1)
% b.FaceColor='b';
% b.EdgeColor='b';
% axis xy tight
% xlim([0 14]); 
% %ylim([0 4])
% hYLabel1 = ylabel('Freq. [x 10^{6}]')
% set(gca,...
%     'FontName'  , 'Helvetica',...
%     'FontSize'  , thisfontsize         ,...
%     'Box'       , 'off'      ,...
%     'OuterPosition', [0.012 0.72 0.515 0.22], ...
%     'TickDir'   , 'out'      ,...
%     'XMinorTick', 'off'      ,...
%     'XTick'     , []         ,...
%     'XTickLabel', []         ,...
%     'Layer'     , 'top'      ,...
%     'YMinorTick', 'on'       ,...
%     'XColor'    , [.0 .0 .0] ,...   
%     'YColor'    , [.0 .0 .0] ,...
%     'LineWidth' , 1         );
% set(hYLabel1,...
%     'FontName', 'Arial',...
%     'FontSize', thisfontsize    );
% 
% subplot(3,3,[6,9])
% b=barh(yedges.*1e5,yn.*1e-6,1), 
% b.FaceColor='r';
% b.EdgeColor='r';
% axis xy tight
% ylim([0.1 2.2]);
% %xlim([0 20]);
% hXLabel1=xlabel('Freq. [x 10^{6}]')
% set(gca,...
%     'FontName'  , 'Helvetica',...
%     'FontSize'  , thisfontsize         ,...
%     'Box'       , 'off'      ,...
%     'OuterPosition', [0.534 0.0010 0.17 0.796],...
%     'TickDir'   , 'out'      ,...
%     'XAxisLocation', 'top'   ,...
%     'XMinorTick', 'off'      ,...
%     'YTick'     , []         ,...
%     'YTickLabel', []         ,...
%     'Layer'     , 'top'      ,...
%     'XMinorTick', 'on'       ,...
%     'XTick'     , 0:20:150      ,...
%     'XColor'    , [.0 .0 .0] ,...   
%     'YColor'    , [.0 .0 .0] ,...
%     'LineWidth' , 1         );
% set(hXLabel1,...
%     'FontName', 'Arial',...
%     'FontSize', thisfontsize    );
% % %xlim([0 10]);
% %hXLabel = xlabel('Absorption index, \beta [x 10^{-7}]');
% %hYLabel = ylabel('Refractive index decrement, \delta [x 10^{-5}] ');
% set(figure(1),'OuterPosition',[402 189 874 720])
% 
% %% Save plots with amplitude
% savedata=1;
% %casefolder=[histogram_path analysis_case '/'];
% %savename=['histogram_VOI_' analysis_case];
% if savedata == 1
%     if ~exist('casefolder','dir'); mkdir(casefolder); end
%     print('-f11','-depsc2', [ casefolder savename '_3D_orientation_all_beta.eps']);
%     print('-f11','-dtiff', [ casefolder savename '_3D_orientation_all_beta.tif']);
%     print('-f8','-depsc2', [ casefolder savename '_3D_orientation_beta.eps']);
%     print('-f8','-dtiff', [ casefolder savename '_3D_orientation_beta.tif']);
%     print('-f9','-depsc2', [ casefolder savename '_histogram_beta.eps']);
%     print('-f9','-dtiff', [ casefolder savename '_histogram_beta.tif']);
%     print('-f10','-depsc2', [ casefolder savename '_bivariate_hist.eps']);
%     print('-f10','-dtiff', [ casefolder savename '_bivariate_hist.tif']);
%     print('-f12','-depsc2', [ casefolder savename '_slice_y_beta.eps']);
%     print('-f12','-dtiff', [ casefolder savename '_slice_y_beta.tif']);
%     print('-f13','-depsc2', [ casefolder savename '_slice_x_beta.eps']);
%     print('-f13','-dtiff', [ casefolder savename '_slice_x_beta.tif']);
%     print('-f14','-depsc2', [ casefolder savename '_slice_z_beta.eps']);
%     print('-f14','-dtiff', [ casefolder savename '_slice_z_beta.tif']);
% end