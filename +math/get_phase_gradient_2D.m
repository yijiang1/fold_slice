% GET_PHASE_GRADIENT_2D get 2D gradient of phase of image IMG. 
% Accept either complex image or just phase 
%
% [d_X, d_Y] = get_phase_gradient_2D(img, step=0, padding=8)
%
% Inputs
%     **img     - stack of complex valued input images 
% *optional*
%     **step    - step used to calculate the central difference, default=0 (analytic expression)
%     **padding - padding around edges to prevent periodic boundary artefacts 
%
% *returns*
%     ++d_X,d_Y - horizontal and vertical phase gradient arrays 
    

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

function [d_X, d_Y] = get_phase_gradient_2D(img, step, padding)

    import utils.*
    import math.*

    if isreal(img)
        img = exp(1i*img);
    end
    if nargin < 2
        step = 0;  % step of the difference (too small will amplify noise)
    end
    if nargin < 3
         padding = 8; % pad arrays to avoid edge issues 
    end
    
    
    % suppress edge issues if phase ramp is not subtracted / there is no
    % air around sample 
    
    if padding > 0
        img = padarray(img,[padding, padding],'symmetric','both'); 
        img = smooth_edges(img, padding, [1,2]);
    end
    
    if step == 0
        % analytic formula (sensitive to noise) but faster 
        % img = img ./ (abs(img) + eps); 
        [d_X, d_Y] = get_img_grad(img);  % img is assumed to be complex 
        d_X = imag(conj(img).*d_X);
        d_Y = imag(conj(img).*d_Y);
    else
        % finite difference based method 
        d_X = angle( imshift_fft_ax(img,-step,2) .* conj( imshift_fft_ax(img,step,2)))/(2*step);
        d_Y = angle( imshift_fft_ax(img,-step,1) .* conj( imshift_fft_ax(img,step,1)))/(2*step);
    end
    if padding > 0
        % remove padding 
        ind = {padding:size(d_X,1)-padding-1,padding:size(d_X,2)-padding-1, ':'};
        d_X = d_X(ind{:}); 
        d_Y = d_Y(ind{:}); 
    end
        

end