close all
addpath(fullfile(pwd,'utils'))
%%
figure1 = figure('Color',[1 1 1],'OuterPosition',[100 100 1800 600]);
%set(gcf,'Outerposition',[100 100 1800 600])
clf()
%set(gca,'color','white')
scale = p.dx_spec; 
%object_roi = object(p.object_ROI{:},:);
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

%pos = pos + Np_o([2,1])/2;
%pos_scales = (pos-Np_o([2,1])/2) .* scale([2,1]); 
pos_scales = pos .* scale([2,1]); 
pos_scales_0 = pos_0 .* scale([2,1]); 
%mean_err = mean(std(pos_err)); 
%range = max(pos) - min(pos); 
%up = 0.02 * min(range) / mean_err; 
%rounding_order = 10^floor(log10(up)); 
%up = ceil(up / rounding_order)*rounding_order;
scatter( pos_scales(:,1),  pos_scales(:,2),'.r')
scatter( pos_scales_0(:,1),  pos_scales_0(:,2),10, '.','MarkerEdgeColor','b')

%hold off 
axis equal xy tight 
range =  [min(pos_scales(:,1)), max(pos_scales(:,1)), min(pos_scales(:,2)), max(pos_scales(:,2))];
axis(range)
title(['\fontsize{16}{\color{blue}initial positions, '...
'\color{red}refined positions}'])
%title('Initial positions: blue. Refined positions: red. Current probe position: Black')
%ylabel('Position [\mum]') %for x-ray
%xlabel('Position [\mum]')
ylabel('Position [A]') %for electron
xlabel('Position [A]')

%%
close all
figure
subplot(2,2,1)
plot(outputs.relative_pixel_scale,'.','MarkerSize',4)
axis tight 
ylabel('Relative pixel scaling correction [-]')
xlabel('Iteration')
title('Scales')
hold off 
grid on 

subplot(2,2,2)
plot(outputs.asymmetry*100,'.','MarkerSize',4)
axis tight 
ylabel('Asymmetry [%]')
xlabel('Iteration')
title('Asymmetry')
grid on 

subplot(2,2,3)
plot(outputs.rotation,'.','MarkerSize',4)
axis tight 
ylabel('Rotation  [deg]')
xlabel('Iteration')
title('Rotation')
grid on 

subplot(2,2,4)
plot(outputs.shear,'.','MarkerSize',4)
axis tight 
ylabel('Shear  [deg]')
xlabel('Iteration')
title('Shear')
grid on 