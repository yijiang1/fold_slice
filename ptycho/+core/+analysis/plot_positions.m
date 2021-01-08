%  PLOT_POSITIONS 
% plot positions of the illumination with respect to the object 
%
% ** p          p structure
%
%
%
%

function plot_positions(p)

    %modified by YJ for electron pty
    if isfield(p,'beam_source') && strcmp(p.beam_source, 'electron')
        unitFactor = 0.1;
        unitLabel = 'nm';
    else %X-ray 
        unitFactor = 1e6;
        unitLabel = '\mum';
    end

    % NOTE: the positions are appended to fig 4 (cf. plot_error_metric)
    numscans = length(p.scan_number);
    scanfirstindex = [1 cumsum(p.numpts)+1]; % First index for scan number
    for ii = 1:numscans
        p.scanidxs{ii} = scanfirstindex(ii):(scanfirstindex(ii+1)-1);
    end

    subplot(2,2,[3 4]); 
    for ii = 1:numscans
        idx = min(ii,p.numobjs);
        positions_centered(p.scanidxs{ii},:) = p.positions(p.scanidxs{ii},:) - p.object_size(idx,:)/2 + p.asize/2;
    end
    cla()
    if p.plot.realaxes
        scale = p.dx_spec*unitFactor;
    else
        % plot positions in pixels, useful for grazing incidence ptycho 
        scale = [1,1];
    end
    hold all 
    for ii = 1:numscans
        plot(positions_centered(p.scanidxs{ii},2).*scale(2), ...
             positions_centered(p.scanidxs{ii},1).*scale(1),...
            'x:','markersize',4); 
    end

    grid on;
    if numscans==1
        title('positions','interpreter','none');
    else
        title('positions (red 1st, blue 2nd)','interpreter','none');
    end
    if p.plot.realaxes
        xlabel(unitLabel);
        ylabel(unitLabel);
    else
        xlabel('pixels');
        ylabel('pixels');
    end
    axis tight equal;
    

end

    
    
