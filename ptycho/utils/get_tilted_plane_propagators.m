% PROP_TILTED_PLANE  Near field propagation into a surface tilted with respect to the beam
%
%    [u_1, H, h_tilted] = prop_free_tilted_plane(u_0, z, lambda, pixel_size, ax=1)
%     returns the propagated wavefield
%    Inputs: 
%       **u_0           stack of images 
%       **rotation      [alpha, beta] - along first, second axis [deg]
%       **lambda        wavelenght [m]
%       **pixel_size    pixel size [m] - in the rotated coordinates, ie. pixel size can be anisotropic 
%    *returns*
%       ++u_1           propagated stack of images 
%
% see utils.prop_free_tilted_plane for more details
    
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

function [fwd_propag, back_propag] = get_tilted_plane_propagators(img_sample, rotation, lambda, pixel_size)
    %%  nearfield propagator into tilted plane 
    import utils.*
    
    fwd_propag = @(x)x; 
    back_propag = @(x)x; 
    
    if all(rotation == 0)
        return
    end
    
    Npix = size(img_sample); 
    
    if any(rotation(1:2)~= 0 )
        % provide tilt of the sample 
        assert(sum(rotation(1:2)~=0) < 2, 'Only rotation along one axis is supported')

        ax = find(rotation(1:2)~=0);

        % propagation distance for each row / column to reach the tilted plane
        % extend * cosd(alpha) * tand(alpha)
        grid = ((-Npix(ax)/2+1):Npix(ax)/2)*pixel_size(min(ax,end)) * sind(rotation(ax));

        [~, H, h_tilted] = prop_free_nonparallel_plane(zeros(Npix(1:2), 'like', img_sample), grid, lambda, pixel_size, ax);

        % precalculate conjuged and transposed matrices 
        H_t = H.'*1;   % enforce copy 
        Hc = conj(H); 
        Hc_t = Hc.'*1; 
        h_tilted_t = h_tilted'*1;   % enforce copy 


        %%%%%%%%%%%%%% propagate the image, see utils.prop_free_tilted_plane for more details  %%%%%%%%%%%%%%%%%%%%%
        if ax == 1
            fwd_propag = @(x)(utils.mtimes_stack(h_tilted,ifft(H_t .* fft(x,[],2),[],2)));
            back_propag =  @(x)(ifft(Hc_t.*fft(utils.mtimes_stack(h_tilted_t,x),[],2),[],2));
        else
            fwd_propag = @(x)(utils.mtimes_stack(ifft(H .* fft(x,[],1),[],1),h_tilted));
            back_propag =  @(x)(ifft(Hc.*fft(utils.mtimes_stack(x,h_tilted_t),[],1),[],1));
        end
    end
    if rotation(3)~= 0
        % rotate image around beam axis 
        fwd_propag =   @(x)fwd_propag(utils.imrotate_ax_fft(x,  rotation(3), 3)); 
        back_propag = @(x)back_propag(utils.imrotate_ax_fft(x, -rotation(3), 3)); 
    end
end





