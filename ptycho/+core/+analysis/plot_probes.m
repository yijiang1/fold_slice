%PLOT_PROBES plot reconstructed probes 
% ** p              p structure
% ** use_display    if false, dont plot results on screen 
%
% *returns*
%  fig - image handle 

function fig3 = plot_probes(p, use_display)
    import utils.rmphaseramp
    import plotting.c2image

    %modified by YJ for electron pty
    if isfield(p,'beam_source') && strcmp(p.beam_source, 'electron')
        unitFactor = 0.1;
        unitLabel = 'nm';
    else %X-ray 
        unitFactor = 1e6;
        unitLabel = '\mum';
    end

    count_plotprb = 1;
    for prmode = 1:p.probe_modes
        for prnum = 1:p.numprobs
        aux = p.probes(:,:,prnum,:);
        E = sum(abs(aux(:)).^2);
            if ~use_display && count_plotprb == 1
                fig3 = plotting.smart_figure('Visible', 'off');
            else
                if count_plotprb == 1 
                     if p.plot.windowautopos && ~ishandle(3) % position it only if the window does not exist
                         fig3 = plotting.smart_figure(3);
                        set(gcf,'Outerposition',[ceil(p.plot.scrsz(4)/p.plot.horz_fact) 1 ceil(p.plot.scrsz(4)/p.plot.horz_fact) ceil(p.plot.scrsz(4)/2)])    %[left, bottom, width, height
                     else
                         fig3 = plotting.smart_figure(3);
                     end
                     clf;
                else
                    set(groot,'CurrentFigure',fig3);
                end
            end
            probe = p.probes(:,:,prnum,prmode); 
            if p.plot.remove_phase_ramp
                probe = rmphaseramp(rmphaseramp(probe,'abs'),'abs');
            end

            subplot(p.plot.subplwinprob(1),p.plot.subplwinprob(2),count_plotprb)
            if ~p.plot.realaxes
                imagesc(c2image(probe));
            else
                iaxis{1} = ([1 p.asize(2)]-floor(p.asize(2)/2)+1)*p.dx_spec(2)*unitFactor; 
                iaxis{2} = ([1 p.asize(1)]-floor(p.asize(1)/2)+1)*p.dx_spec(1)*unitFactor; 
                imagesc(iaxis{:},c2image(probe));
                xlabel(unitLabel)
                ylabel(unitLabel)
            end
            if p.share_probe
                titlestring = sprintf('probe: %s %s', p.plot.prtitlestring, p.plot.extratitlestring);
            else
                titlestring = sprintf('probe: %s %s',p.scan_str{prnum}, p.plot.extratitlestring);
            end
            if p.probe_modes > 1
                Ethis = sum(sum(abs(p.probes(:,:,prnum,prmode)).^2));
                Ethis = Ethis/E;
                titlestring = [titlestring sprintf(' %.1f%%',Ethis*100)];
            end
            title(titlestring,'interpreter','none');
            axis image xy tight

            count_plotprb = count_plotprb + 1;
        end
    end

end
