addpath(strcat(pwd,'/utils/'))

%% 1. Load a .mat reconstruction file w. multiple probe modes 

%% 2. Calculate power percentage of each probe mode
N_probes = size(probe,3);
power = zeros(N_probes,1);
for i = 1:N_probes
    power(i) = mean2(abs(probe(:,:,i,1)).^2);
end
power = power / sum(power);

%% 3. Show probe magnitude and phase
figure
Np = size(probe);
grids = {(-ceil(Np(2)/2):ceil(Np(2)/2)-1)*p.dx_spec(2), ...
     (-ceil(Np(1)/2):ceil(Np(1)/2)-1)*p.dx_spec(1)}; 

for i = 1:N_probes
    mode = mean(mean(probe(:,:,i,1),3),4);
    ax(2*i-1)=subplot(2,N_probes,i);
    RANGE = sp_quantile(abs(mode),[1e-3,1-5e-3], 4)';
    RANGE(2) = max(RANGE(2), RANGE(1)+1e-6);
    amode = abs(mode);
    imagesc3D(grids{:},amode)
    if diff(RANGE)>0;caxis(RANGE); end
    title(sprintf('Mode %i, P:%3.2f', i, power(i)*100))
    %ylabel(sprintf('Amplitude - <%3.2g ; %3.2g>', RANGE))
    axis image xy
    colormap bone 
    set(gca,'TickLength',[0 0])     
    set(gca,'XTick',[],'YTick',[])

    ax(2*i)=subplot(2,N_probes,N_probes+i);
    arg = -angle(stabilize_phase(mode));
    RANGE_arg = sp_quantile(arg,[1e-3,1-1e-3], 4)';
    imagesc3D(grids{:},arg )
    if diff(RANGE_arg)>0;  caxis((RANGE_arg')); end
    %ylabel(sprintf('Phase - <%3.2g ; %3.2g>', RANGE_arg))
    axis image xy
    colormap bone 
    set(gca,'TickLength',[0 0])
    set(gca,'XTick',[],'YTick',[])
end
linkaxes(ax, 'xy')
