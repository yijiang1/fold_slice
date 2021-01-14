% IMSHIFT_CIRC_AX  will apply integer shift that can be different
%     for each frame along axis AX. The shift is applied with !! periodic
%     boundary !!. 
% 
%  img_out = imshift_circ_ax(img, shift, ax)
% 
% Inputs: 
%   **img       - stack of images 
%   **shift     - horizontal / vertical shift in pixels , N*2 vector 
%   **ax        - axis along which the stacked images will be shifted 
% *returns*: 
%   ++img       - shifted image 


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


function img_out = imshift_circ_ax(img, shift, ax)

    shift = round(shift); 
    
    if all(shift == 0)
        img_out=img;
        return
    end
        
    Npix = size(img);
    
    img_out = img; 
         
    
    ind = {':',':',':'}; % assume max 3 dim 
    ax_0 = 1+mod(ax,ndims(img)); % fixed axis 
    for i = 1:Npix(ax_0)
        ind{ax_0} = i;      
        img_out(ind{:}) = circshift(img(ind{:}), shift(i), ax);
    end
    

end

