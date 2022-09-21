function [] = plot_tomo_consistency_shifts(par, theta, shift, total_shift, selected_angles)
%Plot shifts from tomo-consistency alignment for laminography geometry  
% Inputs:
%    **par                     - parameter structure
%    **theta                   - rotation angles
%    **shift                   - shifts given by tomo-consistency alignment
%    **total_shift             - total shifts 
%    **selected_angles         - indices for selected angles to be highlighted in the plots

    import plotting.*
    if ishandle(52); close(52); end
    plotting.smart_figure(52)
    
    subplot1 = subplot(3,1,1);
    [theta_s, ind_sort] = sort(theta);
    shift_total_temp = total_shift + shift;
    shift_total_temp_s = shift_total_temp(ind_sort,:);
    shift_s = shift(ind_sort,:);

    plot(theta_s, shift_s, '.')
    hold on
    if ~isempty(selected_angles)
        plot(theta_s(selected_angles), shift_s(selected_angles,:), 'o', 'MarkerSize', 10)
    end
    legend({'Horizontal shift', 'Vertical shift'}); axis tight ; grid on; xlim(subplot1,[-270.1 90.1]);
    ylabel('Shift [px]'); xlabel('Sorted angles'); title('Shift after self-consistency alignment')
    set(subplot1,'XTick',-270:30:90,'XTickLabel',{'-270','-240','-210','-180','-150','-120','-90','-60','-30','0','30','60','90'});

    subplot2 = subplot(3,1,2);
    plot(theta(ind_sort), shift_total_temp(ind_sort,:)*par.pixel_size(1)/1e-6, '.')
    hold on
    if ~isempty(selected_angles)
        plot(theta_s(selected_angles), shift_total_temp_s(selected_angles,:)*par.pixel_size(1)/1e-6, 'o', 'MarkerSize', 10)
    end
    legend({ 'Horizontal shift', 'Vertical shift'}); axis tight ; grid on; xlim(subplot2,[-270.1 90.1]);
    ylabel('Shift [um]'); xlabel('Angles'); title('Total shift')
    set(subplot2,'XTick',-270:30:90,'XTickLabel',{'-270','-240','-210','-180','-150','-120','-90','-60','-30','0','30','60','90'});

    subplot3 = subplot(3,1,3);
    plot(par.scanstomo, shift_total_temp*par.pixel_size(1)/1e-6, '.')
    hold on
    if ~isempty(selected_angles)
        selected_scans = find(ismember(theta, theta_s(selected_angles)));
        plot(par.scanstomo(selected_scans), shift_total_temp(selected_scans,:)*par.pixel_size(1)/1e-6, 'o', 'MarkerSize', 10)
    end
    legend({ 'Horizontal shift', 'Vertical shift'}); axis tight ; grid on
    ylabel('Shift [um]'); xlabel('Scan numbers'); title('Total shift')
    %drawnow

end


