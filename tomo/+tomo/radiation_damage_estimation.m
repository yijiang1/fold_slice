%  RADIATION_DAMAGE_ESTIMATION Plot SVD filtered curved of the vertical fluctuations 
%  tomo invariant to shown if there was radiation damage 
%
% radiation_damage_estimation(object, par, varargin)
%
% Inputs:
%     **object      - complex valued projections
%     **par         - tomography paramter structure 
%     **varargin    - see the code 

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


function radiation_damage_estimation(object, par, varargin)

    import plotting.*
    import math.*

    
    parser = inputParser;
    parser.addParameter('N_SVD_modes',  1 , @isnumeric )
    parser.addParameter('smoothing', 15 , @isnumeric )
    parser.addParameter('logscale',  true , @islogical ) % threshold for mask estimation 
    parser.addParameter('invariant',  'derivative' , @isstr ) % threshold for mask estimation 
    parser.addParameter('vert_range',  [] , @isnumeric ) % threshold for mask estimation 

    parser.parse(varargin{:})
    r = parser.Results;
    
    if ~isempty(r.vert_range)
        object = object((max(1,r.vert_range(1)):min(end,r.vert_range(end))),:,:);
    end
    
    switch r.invariant
        case 'phase'
            %% standard vertical mass fluctuation
            phase_diff = tomo.get_phase_gradient(object, 2,0.5);
            phase = -tomo.unwrap2D_fft(phase_diff, 2, par.air_gap);
            invariant = max(0,squeeze(sum(phase,2)));
        case 'derivative'
            %% vertical derivative fluctuation
            phase_diff_vert = math.get_phase_gradient_1D(object, 1,1);
            invariant = squeeze(sum(phase_diff_vert,2));
    end

    % smoothing in the time domain to remove outliers 
    mass_filt = medfilt2(invariant,[1,floor(r.smoothing/2)*2+1], 'symmetric'); 
    
    
    figure(5);
    subplot(3,1,1)
    imagesc(invariant)
    axis tight xy
    colormap(franzmap)
    caxis(sp_quantile(invariant, [0.01, 0.99], 10))
    title(sprintf('Vertical mass derivative S%05d - S%05d',par.scanstomo(1),par.scanstomo(end)))
    xlabel('Projection number')
    grid on 
    subplot(3,1,2)
    Navg_slices= 10; % average over last 10 slices 
    mass_filt_resid = mass_filt - median(mass_filt(:,end-Navg_slices:end),2); 
    imagesc(mass_filt_resid)
    axis tight xy
    if r.logscale
        set(gca, 'xscale', 'log')
    end
    colormap bone
    caxis(sp_quantile(mass_filt_resid, [0.01, 0.99], 10))
    title('Filtered change with respect to median')
    xlabel('Projection number')
    grid on 
    
    subplot(3,1,3)
    [U,S,V] = svd(mass_filt); 
    S = S / sqrt(sum(diag(S.^2)));
    plot( V(:,[2:r.N_SVD_modes+1]) * S(2:r.N_SVD_modes+1,2:r.N_SVD_modes+1) )
    if r.logscale
        set(gca, 'xscale', 'log')
    end
    title(sprintf('PCA decomposition - shows changes in the sample, Power=%2.3g%%', 100*sum(diag(S(2:end,2:end)).^2)))
    axis tight
    xlabel('Aprox projection number')
    grid on 
    
    file_png = fullfile(par.output_folder,'Radiation_damage_curve.png');
    disp(['Saving ' file_png])
    
end