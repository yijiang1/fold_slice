% APPLY_3D_APODIZATION Smoothly apodize tomogram to avoid sharp edges and air affecting
% the FRC analysis 
%
% [tomogram,circulo] = apply_3D_apodization(tomogram, rad_apod, axial_apod, radial_smooth)
%
% Inputs:
%     **tomogram - volume to be apodized 
%     **rad_apod - number of pixels to be zeroed from edge of the tomogram 
%     **axial_apod - roughly number of pixels to be zeroed from top / bottom 
%     **radial_smooth - smoothness of the apodization in pixels, default = Npix/10
%     **layer_dim
% Outputs: 
%     ++tomogram - apodized volume 
%     ++circulo -apodization mask 
% MODIFIED BY YJ TO ALLOW UNEVEN SIZES

function [tomogram,circulo] = apply_3D_apodization(tomogram, rad_apod, axial_apod, radial_smooth )
    import utils.*
    [Npix_y,Npix_x,Nlayers] = size(tomogram);
    Npix = max(Npix_y,Npix_x);
    if nargin < 4 
        radial_smooth = Npix/10;
    end

    if nargin < 3
        axial_apod = [];
    end
    if ~isempty(rad_apod) 
        xt = -Npix/2:Npix/2-1;
        [X,Y] = meshgrid(xt,xt);
        radial_smooth = max(radial_smooth,1); % prevent division by zero 
        circulo= single(1-radtap(X,Y,radial_smooth,round(Npix/2-rad_apod-radial_smooth)));
        if Npix_y~=Npix_x 
            circulo= crop_pad( circulo, [Npix_y,Npix_x]);
        end
        tomogram = bsxfun(@times, tomogram, circulo);
    end
    if ~isempty(axial_apod) && Nlayers > 1
        filters = fract_hanning_pad(Nlayers,Nlayers,max(0,round(Nlayers-2*axial_apod)));
        filters = ifftshift(filters(:,1));
        tomogram = bsxfun(@times,tomogram,reshape(filters,1,1,[]));
    end
    
end