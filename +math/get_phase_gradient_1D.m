% GET_PHASE_GRADIENT_1D get 1D gradient of phase of an image stack. 
% Accept either complex image or just phase 
%
% [d_img] = get_phase_gradient_1D(img, ax=2, step=0)
%
% Inputs
%     **img     - stack of complex valued input images 
% *optional*
%     **ax      - axis of derivative, default = 2
%     **step    - step used to calculate the central difference, default=0 (analytic expression)
%
% *returns*
%     ++d_img - phase gradient array
    



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

function d_img = get_phase_gradient_1D(img, ax, step, shift)

    import utils.*
    import math.*

    if isreal(img)
        img = exp(1i*img);
    end
    if nargin < 2
        ax = 2; 
    end
    if nargin < 3
        step = 0.5;  % step of the difference (too small will amplify noise)
    end
    if nargin < 4
        shift = 0;  % perform shift and gradient calculation in single step 
    end
    assert(step >= 0, 'Difference step has to be > 0')
    
    % suppress edge issues if phase ramp is not subtracted / there is no
    % air around sample 
    
    pad_distance = 8; 
    img = padarray(img,circshift([pad_distance,0,0], ax-1),'symmetric','both'); 
    img = smooth_edges(img, pad_distance, ax);

    if step == 0
        % analytic formula (sensitive to noise) but faster 
        img = img ./ (abs(img) + eps); 
        d_img = get_img_grad(img, ax);  % img is assumed to be complex 
        d_img = imag(conj(img).*d_img);
    else
        d_img = angle( imshift_fft_ax(img,-step+shift,ax) .* conj( imshift_fft_ax(img,step+shift,ax)))/(2*step);
    end
    % remove padding 
    ind = circshift({pad_distance:size(d_img,ax)-pad_distance-1,':', ':'},ax-1);
    d_img = d_img(ind{:}); 
    

end