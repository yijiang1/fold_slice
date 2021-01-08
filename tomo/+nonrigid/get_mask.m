% get_mask - estimate support mask for the provided reconstructed volume 
%
% [mask, W_rec] = get_mask(rec_0, mask_threshold, mask_dilate, show_mask)
%
% Inputs:
%    **rec_0                reconstruction  volume 
%    **mask_threshold       relative threshold with respect to the maximum 
%    **mask_dilate          mask dilatation in pixels 
%    **show_mask            true / false if you want to plot the mask 
% Outputs: 
%    ++mask             binary mask 
%    ++W_rec            importance weights for the reconstructioin 

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2018 by Paul Scherrer Institute (http://www.psi.ch)    |
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


function [mask, W_rec] = get_mask(rec_0, mask_threshold, mask_dilate, show_mask)
    
    if nargin < 4
        show_mask = false; 
    end

    %% get mask 
    
    Nlayers = size(rec_0,3); 
    Npix = size(rec_0,1);
    mask_dilate = ceil(mask_dilate); 
    
    mask = rec_0 > mask_threshold * quantile(rec_0(:), 0.99); 
    mask = convn(single(mask), ones(mask_dilate,mask_dilate,mask_dilate, 'single'), 'same') > 1e-3*mask_dilate^3;
    

    %% get importance weighting for difference regions 
    % importance weighting 
    W_rec = gpuArray(single(tukeywin(Npix, 0.2)  .*  tukeywin(Npix, 0.2)' .* reshape(tukeywin(Nlayers, 0.2)',1,1,[])  ));  % avoid edge issues
    W_rec = W_rec .* mask; 
    W_rec = utils.imgaussfilt3_fft(W_rec,mask_dilate/2);
    W_rec = gather(W_rec);

    if show_mask

        plotting.smart_figure(212)
        plotting.imagesc_tomo(mask)
        suptitle('Estimated mask')
        drawnow 

        %% show mask 
        plotting.smart_figure(46)
        subplot(1,2,1)
        plotting.imagesc3D(rec_0.* ~mask, 'init_frame', Nlayers/2)
        colorbar
        axis off image 
        colormap bone 
        title('Example of residuum after applied mask')
        subplot(1,2,2)
        hist(rec_0(1:100:end), 100)
        axis tight 
        drawnow 
    end

end
