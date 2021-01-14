% AFFINE_DEFORM_FFT apply accurate affine deformation on image
% use only for minor corrections !!
%
% img = affine_deform_fft(img, affine_matrix, shift)
%
% Inputs:  
%     **img             - 2D or stack of 2D images 
%     **affine_matrix   - 2x2xN affine matrix 
%     **shift           - Nx2 vector of shifts to be applied 
% *returns*: 
%     ++img             - deformed image 


%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |f
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


function img = imdeform_affine_fft(img, affine_matrix, shift)
    import utils.*
    if nargin < 3
        shift = [];
    end
    if ~isempty(shift)
        img = imshift_fft(img, shift); 
    end
    if ~isempty(affine_matrix)
        if size(affine_matrix,3)>1
            for i=1:size(affine_matrix,3)
                [scale, asymmetry, rotation, shear]  = math.decompose_affine_matrix(double(gather(affine_matrix(:,:,i))));
                if any(abs(scale(:)-1) > 1e-5)
                    img(:,:,i) = imrescale_frft(img(:,:,i), scale, scale.*asymmetry); 
                end
                if any(abs(shear(:)) > 1e-5)
                    img(:,:,i) = imshear_fft(img(:,:,i),shear,1); 
                end
                if any(abs(rotation(:))> 1e-5)
                    img(:,:,i) = imrotate_ax_fft(img(:,:,i),rotation,3); 
                end
            end
            
        else

            [scale, asymmetry, rotation, shear]  = math.decompose_affine_matrix(double(gather(affine_matrix)));
            if any(abs(scale(:)-1) > 1e-5)
                img = imrescale_frft(img, scale, scale.*asymmetry); 
            end
            if any(abs(shear(:)) > 1e-5)
                img = imshear_fft(img,shear,1); 
            end
            if any(abs(rotation(:))> 1e-5)
                img = imrotate_ax_fft(img,rotation,3); 
            end
        end
    end
end