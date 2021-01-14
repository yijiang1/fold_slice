% INTERPOLATEFT_AX Computes interpolated array using 1D Fourier transform, i.e. dirichlet
% interpolation along single axis. Computes the FT and then adjusts the size by zero padding
% or cropping then it computes the IFT. A real valued input may have
% residual imaginary components, which is given by numerical precision of
% the FT and IFT.
%
% imout = interpolateFT_ax(im,outsize,ax, use_fft)
%
% Inputs:
%   **im      - Input complex array
%   **outsize - Output size of array [N pixels]
%   **ax - index of axis along which interpolation is done 
% *optional*
%   **use_fft - if false, assume that im is already fft-transformed, default = true
%
% Outputs:
%   ++imout   - Output complex image
%
% Example: 
% x = randn(10,20,30); 
% x_int = utils.interpolateFT_ax(x, 10, 3) % downsample to 10 pixels along 3rd axis 


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



function [ imout ] = interpolateFT_ax(im,outsize,ax, use_fft)

Nin = size(im);
if ax > ndims(im)
    Nin(ax) = 1; 
end
if nargin < 4
    use_fft = true;
end

Nout = Nin;
Nout(ax) = outsize;

if use_fft
    imFT = fft(im,[],ax);
else
    imFT = im; 
end

centerin = floor(Nin(ax)/2)+1;
centerout = floor(Nout(ax)/2)+1;

center_diff = centerout - centerin;


grid_in = fftshift(1:Nin(ax)); 
grid_in = grid_in(max(-center_diff+1,1):min(-center_diff+Nout(ax),Nin(ax)));
grid_in = {grid_in,':',':',':'}; 
grid_in = circshift( grid_in, ax-1); 

grid_out = [max(ceil(Nout(ax)/2)+1,Nout(ax) - centerin+2):Nout(ax), ...
            1:min(centerin-1, ceil(Nout(ax)/2))]; 
grid_out = {grid_out,':',':',':'}; 
grid_out = circshift( grid_out, ax-1); 



if Nout(ax) > Nin(ax)
    % perform multiplication to keep average values, 
    % multiply the smaller array to save time 
    imFT = imFT*(Nout(ax)/(Nin(ax)));
end

imout = zeros(Nout,'like',im);
imout(grid_out{:}) = imFT(grid_in{:});


if use_fft
    imout = ifft(imout,[],ax);
end

if Nout(ax) < Nin(ax)
    % perform multiplication to keep average values, 
    % multiply the smaller array to save time 
    imout = imout*(Nout(ax)/(Nin(ax)));
end


end
