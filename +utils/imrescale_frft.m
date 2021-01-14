% IMRESCALE_FRFT subpixel accurate image rescaling based on fractional fourier
% transformation (FRFT)
%
% img = imrescale_frft(img, scale_x, scale_y, scale_z)
%
% Inputs:  
%     **img    2D or stack of 2D images 
%     **scale_x - horizontal scaling factor 
%  *optional*
%     **scale_y - vertical scaling factor, if not provided scale_x is used 
%     **scale_z - 3rd axis scaling factor, if not provided, no scaling is
%     used along 3rd axis
% *returns*: 
%     ++img    2D or stack of 2D images scaled by factors scale_x, (scale_y)

    
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


function [img, win] = imrescale_frft(img, scale_x, scale_y, scale_z)

    isReal = isreal(img);
    win = []; 
    if ~isvector(scale_x) && ~isscalar(scale_x)
        error('Inputs scaling is expected as scalar or vector')
    end
    if nargin < 3 && (size(img,1)==size(img,2))
        % 2d version is faster only for many stacked pictures 
        if scale_x > 1
            win = get_window(img, scale_x, 1) .* get_window(img, scale_x, 2);
            img = img .* win; 
        end
        %size(img)
        img = math.fftshift_2D(ifft2(math.fftshift_2D(FRFT_2D(img,scale_x))));
    else
        if nargin < 3
            scale_y = scale_x; 
        end
        if any(scale_y ~= 1)
            img = math.fftshift_2D(ifft(math.fftshift_2D(FRFT_1D(img,scale_y))));
        end
        if any(scale_x ~= 1)
            img = permute(img,[2,1,3]);
            img = math.fftshift_2D(ifft(math.fftshift_2D(FRFT_1D(img,scale_x))));
            img = permute(img,[2,1,3]);
        end
        if nargin > 3
            if any(scale_z ~= 1)
                img = permute(img,[3,2,1]);
                img = math.fftshift_2D(ifft(math.fftshift_2D(FRFT_1D(img,scale_z))));
                img = permute(img,[3,2,1]);
            end
        end
    end
    
    if isReal
        img = real(img);
    end
end
 
function win = get_window(img, scale, ax)
    % apodize window for img to prevent periodic boundary errors 
    win = ones(ceil(size(img,ax)/scale/2)*2,class(img)); 
    win = utils.crop_pad(win, [size(img,ax),1]);
    win = shiftdim(win, 1-ax);
end

function X=FRFT_1D(X,alpha)
    % 1D fractional fourier transformation 
    % See A. Averbuch, "Fast and Accurate Polar Fourier Transform"

    %% it works as magnification lens Claus, D., & Rodenburg, J. M. (2015). Pixel size adjustment in coherent diffractive imaging within the Rayleigh–Sommerfeld regime
    
    %% test plot(abs(fftshift(ifft((FRFT_1D(x,scale))))))
       
    N = size(X,1);
    grid = fftshift(-N:N-1)';

    preFactor = reshape(exp(1i*pi*grid*alpha(:)'),2*N,1,[]);  % perform shift 
    Factor=     reshape(exp(-1i*pi*grid.^2/N * alpha(:)'),2*N,1,[]);  % propagation / scaling
    X=[X; zeros(size(X), class(X))];  % add oversampling 
    X= bsxfun(@times, X, Factor .* preFactor);
    
    % avoid duplication of XX 
    X=fft(X); 
    X = bsxfun(@times, X,fft(conj(Factor)));
    X=ifft(X);
    
    X=bsxfun(@times, X,reshape(Factor .* preFactor,2*N,1,[]));
    X=X(1:N,:,:);
    %% remove phase offset 
    X = bsxfun(@times, X , reshape(exp(-1i*pi*N*alpha/2),1,1,[]));
end

function X=FRFT_2D(X,alpha)
    % 2D fractional fourier transformation 
    % See A. Averbuch, "Fast and Accurate Polar Fourier Transform"

    %% it maybe works as magification lens Claus, D., & Rodenburg, J. M. (2015). Pixel size adjustment in coherent diffractive imaging within the Rayleigh–Sommerfeld regime

    alpha = reshape(alpha,1,1,[]);
       
    N = size(X,1);
    grid = (fftshift(-N:N-1)') * ones(1, 'like', X);
    
    [Xg,Yg] = meshgrid(grid(1:N), grid(1:N));
    preFactor = exp((1i*pi.*alpha)*(-N/2+(Xg+Yg) - (1/N)*(Xg.^2+Yg.^2)));  % perform shift after FFT 
    
    [Xg,Yg] = meshgrid(grid, grid);
    Factor=exp((1i*pi/N(1))*(Xg.^2+Yg.^2) .* alpha);  % propagation / scaling
    Factor = fft2(Factor); 
    
    X=  X .* preFactor;
    
    if length(size(X))==4 %%added by YJ to present errors when using variable probe
        x_tilde = zeros(2*N, 2*N, size(X,3), size(X,4),  'like', X);
        % upsample the X array
        x_tilde(1:N, 1:N,:,:) = X;

        X=fft2(x_tilde);

        X = X .* Factor;

        X=ifft2( X );

        X=X(1:N,1:N,:,:);
        
    else %%length(size(X))==3
    
        x_tilde = zeros(2*N, 2*N, size(X,3),  'like', X);
        % upsample the X array
        x_tilde(1:N, 1:N,:) = X;

        X=fft2(x_tilde);

        X = X .* Factor;

        X=ifft2( X );

        X=X(1:N,1:N,:);
    end
    
    X=X.* preFactor;

end

