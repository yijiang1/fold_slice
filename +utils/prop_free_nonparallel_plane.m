% PROP_FREE_NONPARALLEL_PLANE  Near field propagation into a surface nonparallel with beam
%
%    [U, H, h_tilted] = prop_free_nonparallel_plane(U, z, lambda, pixel_size, ax=1)
%     returns the propagated wavefield
%    Inputs: 
%       **U           stack of images 
%       **z             propagation distance along horizontal axis (vector)
%       **lambda        wavelenght [m]
%       **pixel_size    pixel size [m]
%       **ax            (optional)  axis along which is the tilted plane 
%    *returns*
%       ++U           propagated stack of images 
%
%
%     Example for tilted plane propagation:
%     Npix = size(img);
%     grid = ((-Npix(2)/2+1):Npix(2)/2)*pixel_size * tand(90-grazing_angle); 
%    [img_propag, H, h_tilted] = prop_free_tilted_plane(img, grid, lambda, pixel_size)
% 
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

function [U, H, h_tilted] = prop_free_nonparallel_plane(U, z, lambda, pixel_size, ax)
    %%  nearfield propagator into tilted plane 
    
    if all(z == 0)
        H = 1;
        h = 1; 
        U = U;
        return
    end
    if nargin  < 5
        ax = 1;
    end

    %% standard 1D ASM for vertical axis 
    Np = size(U);
    if ax == 1 
       Np = Np([2,1]);
       pixel_size = pixel_size(min(end,[2,1])); 
    end 
    
    assert(length(z) == Np(2), 'Shift vector z has wrong length')
    
    type = ones(1,'like',U);
    xgrid = type * (-fix(Np(1)/2):ceil(Np(1)/2)-1);
    k = 2 * pi / lambda(1);
    extent = Np .* pixel_size;
    kx =  2 * pi .*xgrid / extent(1) ;
    H = exp(1i*z(:)' .* (kx'.^2) / (2*k));
        
    H = fftshift(H, 1);  % 1D ASM propagator 
    
    %% convolution ASM for horizontal axis
    % fft transform to get real space convolution kernel
    h = ifftshift(fft(H,[],1),1)/(Np(1));
        
    h = utils.crop_pad(h, [Np(2), Np(2)]);

    % apply shift on the convolution propagator for the axis along the beam 
    h_tilted = utils.imshift_fft_ax(h, (1:Np(2))-Np(2)/2-1,1);
    
    % adjustnemnts to make the propagated image consistent with the original one 
    h_tilted = fliplr(rot90(h_tilted));
    H = fliplr(H); 

    %%%%%%%%%%%%%% propagate the image  %%%%%%%%%%%%%%%%%%%%%
    if ax == 1
        % propagate horizontal axis
        U = ifft(H.' .* fft(U,[],2),[],2);
        % propagate vertical axis !! use matrix multiplication !! in realspace
        % to describe the convolution that is needed for propagation to the
        % tilted axis ->   h_tilted * U
        U = utils.mtimes_stack(h_tilted,U);  % convolution applied by matrix multiplication, make it for entire image stack in parallel 
    else
        % propagate horizontal axis
        U = ifft(H .* fft(U,[],1),[],1);
        % propagate vertical axis !! use matrix multiplication !! in realspace
        % to describe the convolution that is needed for propagation to the
        % tilted axis 
        U = utils.mtimes_stack(U,h_tilted);
         
        
    end
    
end