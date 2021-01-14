% IMSHIFT_FFT will apply shift with subpixel accuracy that can be different for  each frame. 
%
% img = imshift_fft(img, x,y, apply_fft = true, weights = [])
%
%  Inputs:
%       **img - input image stack, can be complex valued  
%       **x, y - shifts in number of pixels 
%       **apply_fft , if false, then images will be assumed to be in fourier space 
%       **weights - 0<W<=1 apply importance weighting to avoid noise in low reliability regions 
%  returns: 
%       ++img - shifted image stack 

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



    
function img = imshift_fft(img, x,y, apply_fft, weights)
    
    
    if nargin  < 3
        y = x(:,2);  % x can be either Nx1  or Nx2 vector 
        x = x(:,1);
    end
    if nargin < 4
        apply_fft = true; % if false, assume that img is already fft transformed 
    end
    if nargin < 5
        weights = [];  % weights prevents amplitifaction of noise in low reliability regions 
    end
    if ~isempty(weights) && ~isscalar(weights)
        eps_ = 1e2*eps(ones(1,'like',img)); 
        weights = max(eps_, weights);  % avoid dividing by zero 
    end
    
    if all(x==0) && all(y==0)
        return
    end
    
    if ~isempty(weights) && ~isscalar(weights) && apply_fft
        img = img .* weights; 
    end
    
    if all(x==0)  % shift only along one axis -> faster 
        img = utils.imshift_fft_ax(img, y,1, apply_fft);
    elseif all(y==0)
        img = utils.imshift_fft_ax(img, x,2, apply_fft);
    else
        %%  2D FFT SHIFTING     
        real_img = isreal(img);
        Np = size(img);


        if apply_fft
             img = math.fft2_partial(img);
        end

        xgrid = ifftshift(-fix(Np(2)/2):ceil(Np(2)/2)-1)/Np(2);
        X = reshape((x(:)*xgrid)',1,Np(2),[]);
        X =  exp((-2i*pi)*X);
        img = bsxfun(@times, img,X);
        ygrid = ifftshift(-fix(Np(1)/2):ceil(Np(1)/2)-1)/Np(1);
        Y = reshape((y(:)*ygrid)',Np(1),1,[]);
        Y =  exp((-2i*pi)*Y);
        img = bsxfun(@times, img,Y);

        if apply_fft
            img = math.ifft2_partial(img);
        end
        if real_img
            img = real(img);
        end
    end
    

    if ~isempty(weights) && ~isscalar(weights) && apply_fft
        weights = utils.imshift_fft(weights, x,y);   %% weights needs to be shifted as well
        img = img ./ weights; 
    end
  
  
end
