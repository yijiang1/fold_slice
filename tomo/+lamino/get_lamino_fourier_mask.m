% fft_mask = get_lamino_fourier_mask( Npix, lamino_angle, keep_on_GPU)
% find the missing cone mask based on provided inputs 
% Inputs:
%    **Npix             - (3x1 int) volume size 
%    **lamino_angle     - (scalar), laminography angle from 0 to 90degrees, 90 == classical tomo, it is used to calculate the missing cone  
%    **keep_on_GPU      - (bool) move the mask to GPU and keep it there 
% *returns*
%   ++fft_mask = mask in the fourier space 
%
% Example: 
%   ifftn(fftn(volume).*fft_mask)


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


function fft_mask = get_lamino_fourier_mask( Npix, lamino_angle, keep_on_GPU)
    if nargin < 3
        keep_on_GPU = false; 
    end

    for i = 1:3
        grid{i} = fftshift(linspace(-1,1,Npix(i)))';
        grid{i} = shiftdim(grid{i},1-i);
        if keep_on_GPU, grid{i} = utils.Garray(grid{i}); end
    end
    fft_mask = get_mask(grid{:}, lamino_angle);
end


function fft_mask = get_mask(xgrid, ygrid, zgrid, lamino_angle)
    fft_mask = ceil(atand( abs(zgrid) ./ sqrt(xgrid.^2+ygrid.^2) )) > lamino_angle; 
end
