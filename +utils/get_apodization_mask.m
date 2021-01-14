% GET_APODIZATION_MASK calculate a 2D circular mask fot smoothing tomogram
%
% [circulo] = get_apodization_mask(tomogram, rad_apod, axial_apod, radial_smooth)
%
% Inputs:
%     **tomogram - volume to be apodized 
%     **rad_apod - number of pixels to be zeroed from edge of the tomogram 
%     **radial_smooth - smoothness of the apodization in pixels, default = Npix/10
%     **layer_dim
% Outputs: 
%     ++circulo -apodization mask 
% Written BY YJ based on apply_3D_apodization.m

function [circulo] = get_apodization_mask(Npix, rad_apod, radial_smooth )
    import utils.*
    Npix_y = Npix(1);
    Npix_x = Npix(2);
        
    Npix = max(Npix_y,Npix_x);
    if nargin < 3
        radial_smooth = Npix/10;
    end
    
    if ~isempty(rad_apod) 
        xt = -Npix/2:Npix/2-1;
        [X,Y] = meshgrid(xt,xt);
        radial_smooth = max(radial_smooth,1); % prevent division by zero 
        circulo= single(1-radtap(X,Y,radial_smooth,round(Npix/2-rad_apod-radial_smooth)));
        %size(X)
        if Npix_y~=Npix_x 
            circulo= crop_pad( circulo, [Npix_y,Npix_x]);
        end
    end
end