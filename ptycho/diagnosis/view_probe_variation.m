%view_probe_variation.m
close all
clear
addpath(strcat(pwd,'/utils/'))

% load a .mat reconstruction file w. probe variation correction 
%% calculate probes at each scan positions
Np_p = [size(probe,1),size(probe,2)];
probes = reshape(probe(:,:,1,:),prod(Np_p),[]);
probes = reshape(probes * outputs.probe_evolution', Np_p(1), Np_p(2), []);

%% generate a movie to show scan positions and probes
disp('Generating movie frames')
figure1 = figure('Color',[1 1 1],'OuterPosition',[100 100 1800 600]);
clear M
object_roi = object(p.object_ROI{:},:);
N_sample = 1;
N_frames = 120;
%N_frames = floor(size(outputs.probe_positions(:,1),1)/N_sample);
%scale = p.dx_spec*1e6; % x-ray 
%unit_label = '\mum'; % x-ray 

scale = p.dx_spec; % electron
unit_label = 'A';
for i=1:N_frames
    clf()
    subplot(1,2,1)
    aobject = angle(object);
    range = sp_quantile(angle(object_roi), [1e-3, 1-1e-3],10);
    aobject = (aobject  - range(1)) / (range(2) - range(1));
    Np_o = size(object);
    grids = {(-ceil(Np_o(2)/2):ceil(Np_o(2)/2)-1)*scale(2), ...
             (-ceil(Np_o(1)/2):ceil(Np_o(1)/2)-1)*scale(1)}; 
    imagesc(grids{:}, aobject, [-2, 1]); % reduce contrast 
    colormap bone 
    axis xy
    hold on 
    pos = outputs.probe_positions;
    pos_0 = outputs.probe_positions_0;

    pos_scales = pos .* scale([2,1]); 
    pos_scales_0 = pos_0 .* scale([2,1]); 
    scatter( pos_scales(:,1),  pos_scales(:,2),'.r')
    scatter( pos_scales_0(:,1),  pos_scales_0(:,2),10, '.','MarkerEdgeColor','b')

    scatter( pos_scales(1+(i-1)*N_sample,1),  pos_scales(1+(i-1)*N_sample,2),50,'ok','filled')
    axis equal xy tight 
    range =  [min(pos_scales(:,1)), max(pos_scales(:,1)), min(pos_scales(:,2)), max(pos_scales(:,2))];
    axis(range)
    title(['\fontsize{12}{\color{blue}initial positions, '...
'\color{red}refined positions, \color{black}current probe position}'])
    ylabel(['Position [', unit_label, ']'])

    subplot(1,2,2)

    imagesc_hsv(probes(:,:,1+(i-1)*N_sample),'stabilize_phase',false)
    axis off
    title('primary probe mode')
    drawnow;
    M(i) = getframe(gcf);
end
disp('Done')

%% save movie
disp('Saving movie...')
movieName = strcat('probes_selected.avi');
v= VideoWriter(movieName);

v.FrameRate= 5;
open(v)
writeVideo(v,M);
close(v)
disp('Done')

%% plot OPR modes (U)
N_vp = size(probe,4)-1;
figure
ax(1)=subplot(2,1+N_vp,1);
imagesc_hsv(probe(:,:,1,1),'stabilize_phase',false)
axis xy off
title('Constant mode')

for ii = 2:N_vp+1
    ax(ii)=subplot(2,1+N_vp,ii);
    imagesc_hsv(probe(:,:,1,ii),'stabilize_phase',false)
    axis xy off 
    title(sprintf('Variable mode %i', ii-1))
end

for ii = 1:N_vp+1
    ax(ii+1+N_vp)=subplot(2,1+N_vp,1+N_vp+ii);
    scatter(outputs.probe_positions(:,1)*p.dx_spec(1),outputs.probe_positions(:,2)*p.dx_spec(1),[],outputs.probe_evolution(:,ii));
    
    %axis xy off 
    title(strcat('Average weights:', num2str(mean(outputs.probe_evolution(:,ii)))))
    %title(strcat(num2str(mean(outputs.probe_evolution(:,ii)))))
end

%% save probes as a tiff stack
Np_p = [size(probe,1),size(probe,2)];
probes = reshape(probe(:,:,1,:),prod(Np_p),[]);
probes = reshape(probes * outputs.probe_evolution', Np_p(1), Np_p(2), []);
resampleFactor = round(size(probes,3)/200); %only save ~100 images
probes_temp = probes(:,:,1:resampleFactor:end);
saveName = strcat(reconDir,'probes_Niter',num2str(Niter),'.tiff');
imwrite(convert_to_rgb(probes_temp(:,:,1)), saveName,'tiff')
for i=2:size(probes_temp,3)
    imwrite(convert_to_rgb(probes_temp(:,:,i)), saveName,'tiff', 'WriteMode','append')
end    
