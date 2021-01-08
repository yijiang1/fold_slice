% GET_FSC_FROM_SUBTOMOS calculate Fourier shell correlation between two subtomograms
%
% [resolution FSC T freq n FSC_stats] = get_FSC_from_subtomos(tomograms, FSC_vertical_range, rad_apod, radial_smooth, axial_apod, SNRt,thickring, par)
%
% Inputs:
%   **tomograms - 2x1 cell containing 2 tomograms to be compared 
%   **FSC_vertical_range    - vector of the selected layers for FSC 
%   **rad_apod       - radial apodization of the tomogram volumes 
%   **radial_smooth  - smoothness range of the radial apodization 
%   **axial_apod     - apodizaton along vertical axis 
%   **SNRt           - signal threshold for FRC resolution 
%   **thickring      - thickness of FSC shells 
%   **par            - tomo parameters structure 
% 
% returns: 
%     ++resolution        [min, max] resolution estimated from FSC curve
%     ++FSC               FSC curve values
%     ++T                 Threshold values
%     ++freq              spatial frequencies 
%     ++stat              stat - structure containing other statistics such as
%                         SSNR, area under FSC curve. average SNR, ....
%     ++fsc_path          path to store the FSC curves 

%        
%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%*-----------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing language this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS matlab package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.

function [resolution FSC T freq n FSC_stats, fsc_path] = get_FSC_from_subtomos(tomograms, FSC_vertical_range, rad_apod,radial_smooth,axial_apod,SNRt,thickring,par)

    Npix = size(tomograms{1});
    utils.verbose(struct('prefix', 'FSC'))
    
    
    for ii = 1:2
        auxtomo{ii} = utils.apply_3D_apodization(tomograms{ii}(:,:,FSC_vertical_range(end:-1:1)), rad_apod, axial_apod);
        if ii == 1
            % ignore empty regions 
            tomo_ROI = get_ROI(any(auxtomo{1} ~= 0,3)); 
        end
        auxtomo{ii} = auxtomo{ii}(tomo_ROI{:},:);
    end

    if ishandle(4); close(4); end  % force closing and reopening on the front
    plotting.smart_figure(4)
    subplot(1,2,1)
    img_tmp = rot90(squeeze(tomograms{1}(:,ceil(Npix(1)/2),end:-1:1)),1); 
    imagesc(img_tmp); 
    caxis(math.sp_quantile(img_tmp, [1e-2, 1-1e-2], 10))
    hold on 
    plotting.hline(FSC_vertical_range(1)+axial_apod/2, 'r')
    plotting.hline(FSC_vertical_range(end)-axial_apod/2, 'r')
    plotting.vline(rad_apod+radial_smooth/2, 'b')
    plotting.vline(Npix(1)-rad_apod-radial_smooth/2, 'b')
    hold off 
    % caxis([4,5.3]*1e-3)
    axis xy equal tight
    colormap bone 
    title('Selected FSC range')
    subplot(1,2,2)
    Npix_aux = size(auxtomo{1}); 
    img_tmp = rot90(squeeze(auxtomo{1}(:,ceil(Npix_aux(1)/2),:)),1); 
    imagesc(img_tmp); 
    caxis(math.sp_quantile(img_tmp, [1e-2, 1-1e-2], 10))
    colormap bone 
    hold on 
    plotting.vline(radial_smooth, 'b')
    plotting.vline(Npix_aux(1)-radial_smooth, 'b')
    hold off 
    axis xy equal tight
    title('Input to FSC')
    if par.windowautopos
        win_size = [1000 600]; 
        screensize = get( groot, 'Screensize' );
        set(gcf,'Outerposition',[150 min(270,screensize(4)-win_size(2)) win_size]);
    end
    drawnow 

    if  ~par.online_tomo && ~debug() && strcmpi(input('Accept FSC region (Y/n)? ','s'),'n') 
        utils.verbose(-1,'Manually adjust FSC_vertical_range / axial_apod / rad_apod')
        return
    elseif ~debug()
        try
            fsc_path = fullfile(par.output_folder, sprintf('FSC_region_S%05d_S%05d_%s_freqscl_%0.2f-%s',...
                par.scanstomo(1),par.scanstomo(end),par.filter_type,par.freq_scale,datetime('today')));
            print('-f4','-deps', [fsc_path, '.eps'])
            print('-f4','-dpng', [fsc_path, '.png'])
            utils.verbose(-1,['Saved preview to ', fsc_path])
        catch
            warning('FSC region plot was not saved because figure 4 is missing')
        end
    end


    utils.verbose(-1,'Fourier shell correlation')
    [resolution FSC T freq n FSC_stats] = utils.fourier_shell_corr_3D_2(auxtomo{:},par,'dispfsc',true,'SNRt',SNRt,'auto_binning',true, 'thickring', thickring,'figure_id',45);
    drawnow 

    fsc_path = fullfile(par.output_folder,sprintf('FSC_curve_S%05d_S%05d_%s_freqscl_%0.2f-%s',par.scanstomo(1),par.scanstomo(end),par.filter_type,par.freq_scale,datetime('today')));

    if ~debug()   % ignore in automatic tests 

        print('-f45','-deps', [fsc_path, '.eps'])
        print('-f45','-dpng', [fsc_path, '.png'])
        utils.verbose(-1,['Saved FSC curve to ', fsc_path])

        if par.online_tomo
            fsc_path_online =  sprintf('%s_FSC_tomo.png',par.online_tomo_path); 
            print('-f45','-dpng',fsc_path_online)
            system(sprintf('convert -trim %s %s', fsc_path_online, fsc_path_online));
        end

    end
    utils.verbose(struct('prefix', 'template'))

end
