% FUNCTION  [u_1, H, h, dH] = near_field_evolution(u_0, z, lambda, extent, use_ASM_only)
% Description: nearfield evolution function, it automatically swithch
% between ASM and Fraunhofer propagation 

function [u_1, H, h, dH] = near_field_evolution(u_0, z, lambda, extent, use_ASM_only)


    H = [];
    h = [];
    u_1 = [];
    dH = [];
        

    if nargin < 5
        use_ASM_only = false; 
    end
    
    extent = extent(:)' .* ones(1,2);
    if z == 0
        H = 1;
        u_1 = u_0;
        return
    end
    if z == inf
        return
    end

    Npix = size(u_0);

    xgrid = (0.5+(-Npix(1)/2:Npix(1)/2-1))/Npix(1);
    ygrid = (0.5+(-Npix(2)/2:Npix(2)/2-1))/Npix(2);

    k = 2 * pi / lambda(1);
    
    % Undesamplling parameter
    F = mean( extent.^2 ./   (lambda(1) .* z .* Npix ));
        
    if abs(F) < 1 && ~use_ASM_only
        % farfield propagation 
        warning('Farfield regime, F/Npix=%g', F )
        Xrange = xgrid*extent(1);
        Yrange = ygrid*extent(2);
        [X,Y] = meshgrid(Xrange, Yrange);
        h =  exp(1i*k*z +1i*k/(2*z) * (X'.^2 + Y'.^2));
       
        % this serves as low pass filter for the far nearfield 
        H = ifftshift(fft2(fftshift(h))); 
        H = H / abs(H(end/2+1, end/2+1)); % renormalize to conserve flux in image 
    else
        % standard ASM 
        kx =  2 * pi .*xgrid / extent(1) * Npix(1);
        ky =  2 * pi .*ygrid / extent(2) * Npix(2);
        [Kx, Ky] = meshgrid(kx, ky);
        
        dH =  ( -1i*(Kx'.^2+Ky'.^2)/(2*k) );  
        
        H = exp( 1i*z*sqrt( k^2 - Kx'.^2-Ky'.^2));  % it make it a bit more sensitive to z distance
        h = [];
    end
    
    u_1 = ifft2( bsxfun(@times, ifftshift(H), fft2(u_0)));
end

