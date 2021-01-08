import io.image_read
import io.load_prepared_data
import io.spec_read

clear

addpath ptycho/utils/


scans=[11745:11763];

range=[];

marksaxis_start = 1;    % First markMarks used for scan
marksaxis_step = 0.1;  % Step of marks, change of marks from one scan to the next
usemarksaxis = false;

use_fluorescence = false;
use_pilatus = true;

int=zeros(1,length(scans));
norm=zeros(1,length(scans));
xaxis=zeros(1,length(scans));
marksaxis = [0:numel(xaxis)-1]*marksaxis_step + marksaxis_start;

S=spec_read('~/Data10/specES1/dat-files/','ScanNr',scans);
spectra_sum = [];
for hh=1:length(scans)
    if use_fluorescence
        filename = sprintf('~/Data10/fx123/gccDppConsoleLinux/data/x123_%05d.dat',scans(hh));
        data = image_read(filename);
        if ~isempty(spectra_sum)
            spectra_sum = spectra_sum + sum(data.data,2);
        else 
            spectra_sum = sum(data.data,2);
        end
        if ~isempty(range)
            int(hh)=sum(sum(data.data(range,:),1),2);
        else
            int(hh)=sum(sum(data.data(:,:),1),2);
        end
    end
    if use_pilatus

        prepdata_filename = sprintf('~/Data10/analysis/S%05d/S%05d_data_400x400.h5',scans(hh),scans(hh));
        display(sprintf('loading %s',prepdata_filename))
        [I, fmask, ~] = load_prepared_data(prepdata_filename, true);
        this = sum(I, 3).*any(fmask,3);
        trans(hh) = sum(this(:));
    end
    eaxis(hh)=S{hh}.mokev;
    norm(hh)=S{hh}.bpm4i;
    %data_spectrum = sum(data.data,2);
end




%

if usemarksaxis
    axisplot = marksaxis;
    xlabeltext = 'Marks';
else
    axisplot = eaxis;
    xlabeltext = 'E [keV]';
end

if use_fluorescence
    figure(1)
    subplot(2,1,1)
    semilogy(spectra_sum)
    title('Sum of all spectra')
    
    subplot(2,1,2)
    plot(axisplot,int./norm)
    title('Integrated fluorescence')
    xlabel(xlabeltext)
    
%     figure(3)
%     subplot(2,1,2)
%     line(eaxis,int./norm)
%     %title('Integrated fluorescence')
%     %xlabel(xlabeltext)
%     xlim([eaxis(1) eaxis(end)])
%     ax1 = gca;
%     ax1_pos = ax1.Position; % position of first axes
%     ax2 = axes('Position',ax1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none');
%     line(marksaxis,int./norm,'Parent',ax2,'Color','k')
%     xlim([marksaxis(1) marksaxis(end)])
end

if use_pilatus
    xanes=log(norm./trans);
    
    figure(2)
    title(sprintf('S%05d to S%05d',scans(1),scans(end)))
    subplot(3,1,1)
    plot(axisplot,xanes,'-bo')
    grid on
    
    subplot(3,1,2)
%     plot(axisplot(1:end-1),diff(xanes),'-bo')
    plot(axisplot,1./xanes,'-bo')
    grid on
    
    subplot(3,1,3)
    plot(axisplot, xanes,'-bo')
    xlabel(xlabeltext)
    ylabel('log(bpm4i/pilatus\_transmission)')
    grid on
    
end
