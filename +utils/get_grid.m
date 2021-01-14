% get_grid returns coordinate system for input shape ish and pixel size px
% 
% Example:
%   [g1,g2] = get_grid(512, 29e-9);  
%   returns an fft-shifted coordinate system of size 512x512 with a pixel size of 29 nm

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

function [gx,gy] = get_grid(ish, px)

    if length(ish) == 1 
        sh = [ish ish];
    elseif length(ish) == 2
        sh = ish;
    else
        error('Input shape has to be 1D or 2D')
    end
    
    if length(px) == 1 
        dx = [px px];
    elseif length(ish) ==2
        dx = px;
    else
        error('Pixel size has to be 1D or 2D')
    end
    
    x = fftshift(-sh(2)/2:floor((sh(2)-1)/2))*dx(2);
    y = fftshift(-sh(1)/2:floor((sh(1)-1)/2))*dx(1);
    
    [gx,gy] = meshgrid(x,y);
    
end
