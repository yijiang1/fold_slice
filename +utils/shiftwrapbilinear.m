% out = shiftwraplinear(in,dy,dx)
% Shifts an image with wrap around and bilinear interpolation

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

function out = shiftwrapbilinear(in,dy,dx)

% in = flipud(im2double(imread('cameraman.tif')));
% in = in(1,:);
% %in = in*0;
% %in(100,50) = 1;
% dy = -12;
% dx = -500;

[ny nx] = size(in);
dyfloor = floor(dy);
dxfloor = floor(dx);

x = [1:nx]+dxfloor;
y = [1:ny]+dyfloor;

% Shift integer step (floor of desired)
y = mod(y-1,ny)+1;
x = mod(x-1,nx)+1;

out = in(y,x);


% Subpixel (bilinear)
taux = dx-dxfloor;
tauy = dy-dyfloor;
if (taux~=0)||(tauy~=0)
    indx    = [1:nx];
    indxp1  = [2:nx 1]; 
    indy    = [1:ny];
    indyp1  = [2:ny 1];
    out =   out(indy,indx)*(1-tauy)*(1-taux) + ...
            out(indyp1,indx)*tauy*(1-taux) + ...
            out(indy,indxp1)*(1-tauy)*taux + ...
            out(indyp1,indxp1)*tauy*taux;
end



% figure(100);
% imagesc(out);
% axis xy equal tight;
% colorbar
% figure(101);
% imagesc(real(shiftpp2(in,dy,dx)));
% axis xy equal tight;
% colormap gray
% colorbar