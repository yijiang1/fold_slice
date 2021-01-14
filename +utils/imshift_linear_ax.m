% IMSHIFT_LINEAR_AX  will apply shift that can be different for 
%     each frame along axis ax   
%     + compared to imshift_fft, it does not have periodic boundary 
%     + it is based on linear interpolation, so it can be run fast on GPU 
%     + integer shift is equivalent to imshift_fft (up to the boundary condition)
%     - it needs for-loop for each frame -> it gets slow on GPU for
%       shifting my small images. In that case imshift_fft can be faster. 
%
%   img = imshift_linear(img, x,y, method)
%
%   Inputs:
%       **img       input image / stack of images 
%       **shift     applied shift or vector of shifts for each frame 
%       **ax        axis along which the shift will be performed 
%       **method    choose interpolation method: nearest, {linear}, cubic , circ
%       **extrap_val filling value for the missing regions after interpolation (default=nan)
% 
%   returns: 
%       ++img        shifted image / stack of images 
%
%  see also: utils.imshift_fast, utils.imshift_fft




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



    
function img_out = imshift_linear_ax(img, shift, ax, method, extrap_val)


    if nargin < 4
        method = 'linear';
    end
    if nargin < 5
        extrap_val = nan; 
    end
        
    if all(shift == 0)
        img_out=img;
        return
    end
        
    Npix = size(img);
    
    img = single(img);
    
    img = shiftdim(img, ax-1);

    img_out = img; 
  
    
    ind = {':',':',':'}; % assume max 3 dim 
    ax_0 = 1+mod(ax,ndims(img)); % fixed axis 
    
    if strcmpi(method, 'circ')
        % apply NN shift with circular condition 
        for i = 1:Npix(ax_0)
            ind{ax_0} = i;      
            img_out(ind{:}) = circshift(img(ind{:}), round(shift(i)), ax);
        end
    else
        
        for ii = 1:Npix(ax_0)
            ind{ax_0} = ii;      
            img_out(ind{:}) = interp1(1:size(img,1), img(ind{:}), -shift(ii)+(1:size(img,1)), method, extrap_val);
        end
    
    end
    
    img_out = shiftdim(img_out, ax-1);
    
  
  
end
