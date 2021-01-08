%PLOT_RAW_DATA Simple plotting routine for masked raw data
%
% ** p          p structure
%
% see also: core.initialize_ptycho
function plot_raw_data(p)
import utils.verbose

magnitude = p.fmag .* p.fmask; 

max_intensity = (max(max(magnitude))/p.renorm).^2;

verbose(1, 'Plotting prepared data.')
kk = 1;
title_list = cell(1,size(p.fmag,3));
for jj=1:p.numscans
    for ii = p.scanidxs{jj}
        title_list{kk} = sprintf('Scan S%0.5d - Point (%d)   maximal intensity:%.4g',p.scan_number(jj), ii, max_intensity(kk)); 
        kk = kk +1 ; 
    end
end

% really enforce popup of this figure, it gets very annoying when running somewhere in background
if ishandle(10)
    close(10);
end
fig = figure(10); 

if ~p.fourier_ptycho
    plt_fnct = @(x)(log10(0.1+math.fftshift_2D(x / p.renorm).^2));
else
    plt_fnct = @(x)((x / p.renorm).^2);
end

plotting.imagesc3D(magnitude, 'title_list', title_list, 'fnct', plt_fnct);
if ~p.fourier_ptycho
    caxis([-1, log10(max(max_intensity))])
else
    caxis(math.sp_quantile((magnitude/p.renorm).^2, [1e-6 1-1e-6], 10))
end

c = colorbar; 
ylabel(c, 'log10 counts')
ax = fig.CurrentAxes;
axis(ax, 'xy', 'equal', 'image')
colormap(ax, 'plotting.franzmap')
if p.plot.windowautopos
    horiz_fact = 2.5; 
    if check_option(p.plot, 'object_spectrum')
        pos = 3;
    else
        pos = 2;
    end
    set(gcf,'Outerposition',[ceil(min(p.plot.scrsz(3)-ceil(p.plot.scrsz(4)/2), ceil(p.plot.scrsz(4)*pos/horiz_fact)))  ceil(p.plot.scrsz(4)/2) ceil(p.plot.scrsz(4)/horiz_fact) ceil(p.plot.scrsz(4)/2)])    %[left, bottom, width, height
end

ax.play(ax);

end

