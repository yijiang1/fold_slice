%  PLOT_PROBES_AT_DETECTOR 
% plot reconstructed probes propagated to the detector 
%
% ** p                  p structure
% ** use_display        if false, do now show plots 
%
% *returns*
%  ++fig - image handle 
%
%

function fig5 = plot_probes_at_detector(p, use_display)

    count_plotprb = 1;
    for prmode = 1:p.probe_modes
        for prnum = 1:p.numprobs
        aux = p.probes(:,:,prnum,:);
        E = sum(abs(aux(:)).^2);
            if ~use_display && count_plotprb == 1
                fig5 = plotting.smart_figure('Visible', 'off');
            else
                if count_plotprb == 1
                     if p.plot.windowautopos && ~ishandle(5) % position it only if the window does not exist
                         fig5 = plotting.smart_figure(5);
                         set(gcf,'Outerposition',[ceil(p.plot.scrsz(4)*2/p.plot.horz_fact) 1 ceil(p.plot.scrsz(4)/p.plot.horz_fact) ceil(p.plot.scrsz(4)/2)])    %[left, bottom, width, height
                     else
                         fig5 = plotting.smart_figure(5);
                     end
                     clf;
                else
                    set(groot,'CurrentFigure',fig5);
                end
            end
            subplot(p.plot.subplwinprob(1),p.plot.subplwinprob(2),count_plotprb)
            af_probe = abs(fftshift(fft2(p.probes(:,:,prnum,prmode)))).^2; 
            max_af_probe = max(af_probe(:)); 
            if isfield(p, 'renorm')
                af_probe = af_probe / single(p.renorm).^2; 
            end
            if ~p.plot.realaxes
                imagesc(log10(1e-2*max_af_probe+af_probe));
            else
                imagesc(([1 p.asize(2)]-floor(p.asize(2)/2)+1)*p.ds*1e3,([1 p.asize(1)]-floor(p.asize(1)/2)+1)*p.ds*1e3,log10(1e-2*max_af_probe+af_probe));
                xlabel('mm')
                ylabel('mm')
            end
            if p.share_probe
                titlestring = sprintf('log10 FFT probe: %s %s', p.plot.prtitlestring, p.plot.extratitlestring);
            else
                titlestring = sprintf('log10 FFT probe: %s %s',p.scan_str{prnum}, p.plot.extratitlestring);
            end
            if p.probe_modes > 1
                Ethis = sum(sum(abs(p.probes(:,:,prnum,prmode)).^2));
                Ethis = Ethis/E;
                titlestring = [titlestring sprintf(' %.1f%%',Ethis*100)];
            end
            title(titlestring,'interpreter','none');
            axis image xy tight
            colormap(plotting.franzmap)
            colorbar
            count_plotprb = count_plotprb + 1;
        end
    end
    
end

