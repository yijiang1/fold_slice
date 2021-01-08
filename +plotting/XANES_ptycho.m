%% Loading data
addpath ..
addpath ../ptycho/utils/
addpath ../ptycho/

clear

% close all;

scans =  [3917:3997]; % scan numbers

for ii=1:numel(scans)
    ptycho_filename = find_ptycho_filename('../../',scans(ii),[],'recons');
    if ~iscell(ptycho_filename)
        filename_list{ii}=utils.abspath(ptycho_filename);
    else
        filename_list{ii}=utils.abspath(ptycho_filename{end});
    end
end


    p0=io.HDF.hdf5_load(filename_list{1},'/reconstruction/p');
    
    Nthreads = 6;
    
    dims_ob=p0.object_size;
    pixsize=p0.dx_spec(1);
    size_crop=round((p0.object_size-p0.asize)/2)*2;

    
    proj_file_names  = reshape(filename_list, 1,[]); 
    
    object_block = io.ptycho_read(Nthreads, 'single', dims_ob, '/reconstruction/object', proj_file_names); 
    object_block=permute(object_block,[2 1 3]);
    size_now=size(object_block);
    pixsize=zeros(numel(scans),1);
    
    for ii=1:numel(scans)
        energies(ii)=io.HDF.hdf5_load(filename_list{ii},'/reconstruction/p/energy');
        pixsize(ii,:)=io.HDF.hdf5_load(filename_list{ii},'/reconstruction/p/dx_spec');
    end
    pixsize = pixsize(:,1);
    ROI = {round(size_now(1)/2) + (-size_crop(1)/2:size_crop(1)/2),round(size_now(2)/2)+(-size_crop(2)/2:size_crop(2)/2)};
    object_crop=utils.imrescale_fft(object_block(ROI{:},:),pixsize(ii)/max(pixsize));

	    
    amp=abs(object_crop);
    phase=math.unwrap2D_fft2(object_crop,10,0,[],-1);


%% Selecting the nonair part
roi_chosen = {1:379,70:600}; 
% If the sample is in the middle of FOV, include the entire object so only
% air is outside this ROI.
% If air region is in the middle of FOV, include only air.

obj_or_air=1; % =1 if object is in middle, =0 if air is in middle.

figure(5); clf;
plotting.imagesc3D(phase); 
colormap bone; axis xy equal tight;
caxis([-2*pi 0.2]);
hold on;
plotting.vline(roi_chosen{2}(1),'-b');
plotting.vline(roi_chosen{2}(end),'-b');
plotting.hline(roi_chosen{1}(1),'-b');
plotting.hline(roi_chosen{1}(end),'-b');

drawnow;

figure(6); clf;
plotting.imagesc3D(amp); 
colormap bone; axis xy equal tight;
caxis(math.sp_quantile(amp,[0.01 0.99],10));
hold on;

plotting.vline(roi_chosen{2}(1),'-b');
plotting.vline(roi_chosen{2}(end),'-b');
plotting.hline(roi_chosen{1}(1),'-b');
plotting.hline(roi_chosen{1}(end),'-b');

size_t=size(amp);

if obj_or_air==1
    mask=ones(size_t(1),size_t(2));
    mask(roi_chosen{:})=0;
else
    mask=zeros(size_t(1),size_t(2));
    mask(roi_chosen{:})=1;
end

for ii=1:numel(scans)
    trans = abs(object_crop(:,:,ii));
    reg(ii)=mean(mean(trans(mask==1)));
    trans = trans/reg(ii);
    transmission (:,:,ii) = trans;
end

%% Alignment of projections
margin = min(roi_chosen{2}(1),size(phase,2)-roi_chosen{2}(2));
phase_unwrap=math.unwrap2D_fft2(object_crop,margin,0,[],-1);

shift = utils.find_shift_fast_2D(phase_unwrap, ref, 0.1); 
phase_unwrap_aligned = utils.imshift_fft(phase_unwrap, -shift); 
transmission_aligned = utils.imshift_fft(transmission, -shift); 

%ref=mean(phase_unwrap,3);
%for ii=1:numel(scans)
%   curr=phase_unwrap(:,:,ii);
%    [s1,s2,delta]= utils.registersubimages_2(curr,ref,[],[],[],[],100,0,1);
%    translate(:,ii)=delta(:);
%    phase_unwrap(:,:,ii)=utils.shiftwrapbilinear(curr,-delta(1),-delta(2));
%    transmission_aligned(:,:,ii)=utils.shiftwrapbilinear(transmission(:,:,ii),-delta(1),-delta(2));   
%end

figure(5); clf; 
plotting.imagesc3D(phase_unwrap_aligned); 
colormap bone; axis xy equal tight;
caxis([-2*pi 0.2]);
hold on;
plotting.vline(roi_chosen{2}(1),'-b');
plotting.vline(roi_chosen{2}(end),'-b');
plotting.hline(roi_chosen{1}(1),'-b');
plotting.hline(roi_chosen{1}(end),'-b');
drawnow;

figure(6); clf; 

plotting.imagesc3D(transmission_aligned); 
colormap bone; axis xy equal tight;
caxis(math.sp_quantile(transmission_aligned,[0.01 0.99],10));
hold on;
plotting.vline(roi_chosen{2}(1),'-b');
plotting.vline(roi_chosen{2}(end),'-b');
plotting.hline(roi_chosen{1}(1),'-b');
plotting.hline(roi_chosen{1}(end),'-b');

%% XANES analysis

roi = {50:250;200:400}; %{15:28;22:35}

figure(5); clf; 
plotting.imagesc3D(phase_unwrap_aligned); 
colormap bone; axis xy equal tight;
caxis([-2*pi 0.2]);
hold on;
plotting.vline(roi{2}(1),'-r');
plotting.vline(roi{2}(end),'-r');
plotting.hline(roi{1}(1),'-r');
plotting.hline(roi{1}(end),'-r');

plotting.vline(roi_chosen{2}(1),'-b');
plotting.vline(roi_chosen{2}(end),'-b');
plotting.hline(roi_chosen{1}(1),'-b');
plotting.hline(roi_chosen{1}(end),'-b');

figure(6); clf; 
plotting.imagesc3D(transmission_aligned); 
colormap bone; axis xy equal tight;
caxis(math.sp_quantile(transmission_aligned,[0.01 0.99],10));
hold on;
plotting.vline(roi{2}(1),'-r');
plotting.vline(roi{2}(end),'-r');
plotting.hline(roi{1}(1),'-r');
plotting.hline(roi{1}(end),'-r');

plotting.vline(roi_chosen{2}(1),'-b');
plotting.vline(roi_chosen{2}(end),'-b');
plotting.hline(roi_chosen{1}(1),'-b');
plotting.hline(roi_chosen{1}(end),'-b');


if ~isempty(roi)
    xanes_phase = squeeze(mean(mean(phase_unwrap_aligned(roi{:},:))));
    xanes_amp   = squeeze(mean(mean(transmission_aligned(roi{:},:))));
else
    xanes_phase = squeeze(mean(mean(phase_unwrap_aligned)));
    xanes_amp   = squeeze(mean(mean(transmission_aligned)));
end
figure(7); plot(energies,-xanes_phase(:),'ro-'); %.*energies(:).^2);
title('phase')
figure(8); plot(energies,xanes_amp(:),'ro-');
title('amp')
figure(9); plot(energies(1:end-1),diff(xanes_amp(:)));
title('amp\_diff')

%%
% estep=[3:15];
% figure(7); hold off;
% plot(energies,xanes_phase(:)); hold on;
% plot(energies(estep),xanes_phase(estep),'ro');
% figure(8); hold off;
% plot(energies,xanes_amp(:)); hold on;
% plot(energies(estep),xanes_amp(estep),'ro');
