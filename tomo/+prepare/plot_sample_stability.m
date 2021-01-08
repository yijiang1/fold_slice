%    PLOT_SAMPLE_STABILITY   plot time evolution of the sample stability estimated from the
%       OMNY intererometers, useful to detect unexpected behaviour / setup vibrations 
%
%   poor_projections = plot_sample_stability(par, scanstomo,plot_stability, vibrations_threshold)
%
% Inputs:
%   ++scanstomo           - number of the scans that will be checked 
%   ++plot_stability      - (bool, default = true) 
%   ++vibrations_threshold - (default = pixel_size), acceptable level of vibrations in nanometers 
% *returns* 
%   **poor_projections    - (bool array), true if the projection stability is worse than vibrations_threshold


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



function poor_projections = plot_sample_stability(par, scanstomo, plot_stability, vibrations_threshold)

    if nargin < 3
        plot_stability = true; 
    end
    if nargin < 4
        vibrations_threshold = par.pixel_scale; 
    end
    
    num_proj = length(scanstomo);
    if num_proj == 0
        poor_projections = []; 
        return
    end
    
    for ii = 1:num_proj
        utils.progressbar(ii,num_proj)
        out = beamline.read_omny_pos(sprintf(par.omnyposfile, scanstomo(ii)));
        std_err(ii,1) = quantile(out.Stdev_x_st_fzp,0.8);
        std_err(ii,2) = quantile(out.Stdev_y_st_fzp,0.8);
        d = dir(sprintf(par.omnyposfile, scanstomo(ii)));
        scan_time(ii) = d.datenum;
    end

    if plot_stability
        plotting.smart_figure(5457)

        % Create the first axes
        hax1 = axes();

        % Plot something here
        line(scan_time,std_err*1e3,'color', 'white');
        datetick('x','HH:MM')
        xlim([min(scan_time), max(scan_time)])
        grid on
        ylabel('Average vibrations [nm]')
        xlabel('Scan time')

        % Create a transparent axes on top of the first one with it's xaxis on top
        % and no ytick marks (or labels)
        hax2 = axes('Position', get(hax1, 'Position'), ...  % Copy position
                    'XAxisLocation', 'top', ...             % Put the x axis on top
                    'YAxisLocation', 'right', ...           % Doesn't really matter
                    'xlim', [scanstomo(1),scanstomo(end)], ...                     % Set XLims to fit our data
                    'Color', 'none', ...                    % Make it transparent
                    'YTick', []);                           % Don't show markers on y axis
        % Plot data with a different x-range here

        hplot2 = line(scanstomo,std_err*1e3, 'Parent', hax2);
        legend({'Horizontal','Vertical'}, 'Location', 'Best')
        xlabel(hax2, 'Scan number')
    
        % Link the y limits and position together
        linkprop([hax1, hax2], {'ylim', 'Position'});
        
        plotting.hline(vibrations_threshold)
        
    end
    
    % estimate poor projections 
    
    poor_projections = any(std_err > 1e6*vibrations_threshold,2); 
    

end
