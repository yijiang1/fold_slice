% show_rigid_corrections - plot recovered shifts for each of the projection
%
%   show_rigid_corrections(rec,  sinogram_shifted, err,shift_all, angles, iter, par)
%
% Inputs:
%    **rec          current reconstruction 
%    **sinogram_shifted   sinogram with already applied position shifts 
%    **err          projection space (sinogram - model) error 
%    **shift_all    reconstructed projections 
%    **angles       angles of the projections 
%    **iter         current iteration 
%    **par          parameter structure of the nonrigid tomo 


%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
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


function show_rigid_corrections(rec,  sinogram_shifted, err,shift_all, angles, iter, par)
        import utils.*
        import math.*

        [Nlayers,~,~] = size(sinogram_shifted);
        
        [angles, ind_sort] = sort(angles); 
        sinogram_shifted = sinogram_shifted(:,:,ind_sort); 
        err = err(:,ind_sort); 
        shift_all = shift_all(:,ind_sort,:); 

        
        verbose(1,'Plotting')
        
        figure(5464)
        clf()
        subplot(2,3,1)

        imagesc(squeeze(sinogram_shifted(ceil(Nlayers/2),:,:))');  
        axis off 
        colormap bone 
        title('Corrected sinogram')
        
        subplot(2,3,2)
        if iter > 1
            hold on 
            plot(angles, (shift_all(iter,:,1)-shift_all(iter-1,:,1))*par.binning, 'r')
            plot(angles, (shift_all(iter,:,2)-shift_all(iter-1,:,2))*par.binning, 'b')
            hold off 
            legend({'horiz', 'vert'})
        end
        title('Current position update')
        xlim([min(angles), max(angles)])
        ylabel('Shift [px]')
        xlabel('Angle [deg]')

        subplot(2,3,3)
        hold on 
        plot(angles,shift_all(iter,:,1)*par.binning, 'r')
        plot(angles,shift_all(iter,:,2)*par.binning, 'b')
        hold off 
        title('Total position update')
        legend({'horiz', 'vert'})
        ylabel('Shift [px]')
        xlim([min(angles), max(angles)])
        xlabel('Angle [deg]')

        
        subplot(2,3,4)
        Nlayers = size(rec,3);
        plotting.imagesc3D(rec, 'init_frame', ceil(Nlayers/2))
        caxis(gather(math.sp_quantile(rec(:,:,ceil(Nlayers/2)), [0.01,0.99], 1)));
        axis off image
        title('Current reconstruction')
        colormap bone 
        subplot(2,3,5)
        hold on
        plot(err)
        plot(mean(err,2), 'k', 'LineWidth', 3);
        hold off 
        grid on 
        axis tight
        xlim([1,iter+1])
        set(gca, 'xscale', 'log')
        set(gca, 'yscale', 'log')
        title('MSE evolution')
        xlabel('Iteration')
        ylabel('Mean square error')
    
          
        subplot(2,3,6)
        hold on 
        plot(angles, err(end,:), 'k.') 
        hold off 
        if any(~par.valid_angles)
            legend({'errors', 'ignored'})
        end
        title('Current error')
        xlim([min(angles), max(angles)])
        xlabel('Angle [deg]')

        drawnow 

        
end
