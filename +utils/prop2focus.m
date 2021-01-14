%PROP2FOCUS propagate img to focus
%   prop2focus uses 'phase detection autofocus' to find the focus
%   
%   img...          complex-valued object
%   lam...          wavelength
%   dx...           pixel size
%
%   optional parameters
%   d_start...      initial guess of propagation distance to focus
%   fov...          crop to fov and apodize edges

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

function [ d, img_prop ] = prop2focus(img, lam, dx, varargin)
import utils.*

img_sz = size(img);

% defaults
fov = img_sz(1)*0.6;
d_start = 0;

% parse the variable input arguments vararg = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch lower(name)
        case 'd_start'
            d_start = value;
        case 'fov'
            fov = value;
        otherwise
            error('Unknown parameter %s', name)
    end
end

% prepare fov mask with apodization
fov_mask = fftshift(fract_hanning_pad(img_sz, round(fov*1.2), round(fov)));
img = img.*fov_mask;


% create masks for autofocus
mask = fract_hanning_pad(img_sz,round(img_sz(1)/5),round(img_sz(1)/5*0.9));

mask1 = abs(shiftpp2(fftshift(mask),round(img_sz(1)/20),0));
mask2 = abs(shiftpp2(fftshift(mask),-round(img_sz(1)/20),0));


% minimize difference between the 2 images
fun = @(d)sum(sum(abs(abs(fft2(ifftshift(fftshift(fft2(prop_free_nf(img,lam,d,dx))).*mask2)))-(abs(fft2(ifftshift(fftshift(fft2(prop_free_nf(img,lam,d,dx))).*mask1)))))));

d = fminsearch(fun,d_start);

% propagate to focus
img_prop = prop_free_nf(img, lam, d, dx);


end
