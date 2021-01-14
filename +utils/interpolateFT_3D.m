% INTERPOLATEFT_3D Computes 3D interpolated image using Fourier transform, i.e. dirichlet
% interpolation. Computes the FT and then adjusts the size by zero padding
% or cropping then it computes the IFT. A real valued input may have
% residual imaginary components, which is given by numerical precision of
% the FT and IFT.
%
%  imout = interpolateFT_3D(im,outsize, fourier_mask)
%
%   Inputs
%       **im      - Input real/complex 3D volume
%       **outsize - Output size of array [ny nx nz]
%       **fourier_mask - if provided, apply mask in fourier space. ifftn( mask * fftn(im))
%   *returns*
%       ++imout   - Output real/complex volume

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

function [ imout ] = interpolateFT_3D(im,outsize, fourier_mask)

Nout = outsize;
Nin = size(im);

if all(Nout == Nin) && nargin < 3
    imout = im; 
    return
end



imFT = fftshift(fftn(im));

imout = utils.crop_pad_3D(imFT, outsize);

if nargin > 2
    % if provided, apply fourier mask, NOT FFTSHIFTED !!
    imout = imout .* fourier_mask;    
end

imout = ifftn(ifftshift(imout))*(Nout(1)*Nout(2)*Nout(3)/(Nin(1)*Nin(2)*Nin(3)));

if isreal(im)
    imout = real(imout); 
end





