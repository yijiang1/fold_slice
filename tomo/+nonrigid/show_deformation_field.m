% show_deformation_field - plot reconstructed deformation vector field 
%
%   show_deformation_field(rec_avg, deform_tensors, apodize_radial, Nsvd, binning, upscale_arrows, slice_axis, down_DVF)
%
% Inputs:
%    **rec_avg            optimal reconstuction 
%    **deform_tensors     reconstructed DVF 
%    **apodize_radial     apply radial appodization to crop artefacts around 
%    **Nsvd               number of SVD modes to be plotted 
%    **binning            currenlty used binning (used for scaling)
%    **upscale_arrows     (scalar) scaling constant for the plotted arrows 
%    **slice_axis         axis along which the reconstruction will be sliced and plotted    
%    **down_DVF           (int) donsample DVF to make the arrows more sparse in the plot
%


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


function show_deformation_field(rec_avg, deform_tensors, apodize_radial, Nsvd, binning, upscale_arrows, slice_axis, down_DVF)
    % show deformation vector field 

    Nblocks = length(deform_tensors); 
    
    for kk = 1:3
        for ll = 1:Nblocks+1
            if ll <= Nblocks
                shift_3D_mat(:,:,:,kk,ll) = deform_tensors{ll}{1,kk};
            else
                shift_3D_mat(:,:,:,kk,ll) = deform_tensors{ll-1}{min(end,2),kk};
            end
        end
    end


    [~,mask_small] = utils.apply_3D_apodization(deform_tensors{1}{1}, apodize_radial); 

    

    % apply mask on the results 
    shift_3D_mat = shift_3D_mat .* mask_small; 
    
    % swap dimension to show the right plane 
    switch slice_axis 
        case 1
            rec_avg = rot90(permute(rec_avg, [2,3,1]),1);
            mask_small = rot90(permute(mask_small, [2,3,1]),1);
            shift_3D_mat = rot90(permute(shift_3D_mat, [2,3,1,4,5]),1);
        case 2
            rec_avg = rot90(permute(rec_avg, [1,3,2]),1);
            mask_small = rot90(permute(mask_small, [1,3,2]),1);
            shift_3D_mat = rot90(permute(shift_3D_mat, [2,3,1,4,5]),1);
        case 3
    end
    
    
    [Nx,Ny,Nlayers] = size(rec_avg);
    Nps = size(shift_3D_mat); 
    mesh_2D = {1:down_DVF:Nps(1),1:down_DVF:Nps(2)};
    frame_s = ceil(Nps(3)/2);
    frame = ceil(Nlayers/2); 
    
    % calculate SVD 
    shift_3D_mat = reshape(shift_3D_mat,[],Nblocks+1);
    [U,S,V] = math.fsvd(shift_3D_mat, Nsvd);
    U = reshape(U, [Nps(1:3), 3, Nsvd]);

    U = U * sign(mean(V(:,1))); 
    V = V * sign(mean(V(:,1))); 
    if slice_axis == 3
        mean_amp = sqrt(mean(math.mean2(abs(U).^2 .* mask_small) ./ math.mean2(mask_small),3)); 
    else
        mean_amp = sqrt(mean(math.mean2(abs(U).^2))); 
    end
    mean_amp = squeeze(mean_amp(1,1,1,:,1)); 
    mean_amp = mean_amp .* S(1) * V(3,1); 
    mean_amp = mean_amp .* binning;
    fprintf('Mean deformation x:%3.2gpx y:%3.2gpx z:%3.2gpx \n',mean_amp )
    [~,S_tmp,~] = math.fsvd(shift_3D_mat, min(size(shift_3D_mat,2),Nsvd+10));
    fprintf(['Relative power of the modes:', repmat(' %3.3g%%, ',1,size(S_tmp,1)) , ' \n'], diag(S_tmp ./ sum(S_tmp(:)))*100 )


    figure(545)
    for mode = 1:Nsvd
        ax(mode) = subplot(2,Nsvd,mode);
        for kk = 1:3
            Q{mode,kk} = U(:,:,:,kk,mode) .* S(mode,mode);
            Q{mode,kk} = utils.imgaussfilt3_fft(Q{mode,kk}, down_DVF);
        end
        ygrid = ((1:Nps(1))-0.5)/Nps(1)*Nx; 
        xgrid = ((1:Nps(2))-0.5)/Nps(2)*Ny; 
        [x,y] = meshgrid(xgrid, ygrid); 
        img = rec_avg(:,:,frame); 
        img = min(1,img / math.sp_quantile(rec_avg(:), 0.95,5)); 
        imagesc(1-img, [-1,1]);
        colormap bone 
        hold all

        quiver(x(mesh_2D{:}),y(mesh_2D{:}),Q{mode,2}(mesh_2D{:},frame_s)*upscale_arrows, ...
            Q{mode,1}(mesh_2D{:},frame_s)*upscale_arrows,0,'Linewidth',2); 
        axis off  image 
        hold off 
        title(sprintf('%i. PCA of DVF field\n %ix upscaled',mode, upscale_arrows))
        subplot(2,Nsvd,Nsvd + mode)
        plot(V(:,mode))
        grid on 
        axis tight 
        xlabel('Interpolation node id')
        ylabel('Normalized evolution')
    end
    linkaxes(ax, 'xy')
    plotting.suptitle('Singular value decomposition of the DVF evolution')
    % end

    figure(45545)
    subplot(1,2,1)
    plotting.imagesc3D(rec_avg, 'init_frame', size(rec_avg,3)/2)
    axis off image ; colormap bone 
    colorbar  
    title('Reconstruction example')
    subplot(1,2,2)
    deform = Q{1,1}*binning;
    plotting.imagesc3D(deform, 'init_frame', size(deform,3)/2)
    caxis(gather(math.sp_quantile(deform, [0.001, 0.999],5)))
    axis off image ; colormap bone 
    title('Vertical deformation vector field')
    plotting.suptitle('1th PCA vector, horizontal cut')

    drawnow 

end



