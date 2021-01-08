% UNWRAP2D_FFT  simple and very fast 1D phase unwrapping applied for each
% slice of the provided image stack 
%
% phase = unwrap2D_fft(img, axis, empty_region, step)
%
% Inputs: 
%   **img           - (2D or 3D image stack) either complex valued image or real valued phase gradient 
%   **axis          - (scalar, int) axis along which the gradient is taken 
% *optional*
%   **empty_region  - 2x1 or 1x1 vector, size of empty region assumed around edges for phase offset removal, default=[] 
%   **step          - used to calculate finite difference gradient, 0 = analytical (default) expression
%
% *returns*
%   ++phase         - unwrapped phase with phase ramp removed 
%   ++phase_diff    - if img is complex array, phase difference can be also returned 
%   ++residues      - calculated binary map of phase residuas 
%
% See also: utils.findresidues, utils.remove_sinogram_ramp


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




function [phase, phase_diff, residues] = unwrap2D_fft(phase_diff, axis, empty_region, step)

    import math.*
    import utils.*

    if nargin < 4
        step = 0;
    end
    if nargin < 3
        empty_region = [];
    end
    if nargout > 2 && ~isreal(phase_diff)  % if input is complex array
        residues = abs(findresidues(phase_diff)) > 0.1; 
    else
        residues = []; 
    end

    if ~isreal(phase_diff)
        phase_diff = get_phase_gradient_1D(phase_diff, axis, step);
    end
    phase = real(get_img_int_1D(phase_diff,axis));

    if ~isempty(empty_region) && axis == 2
       phase = remove_sinogram_ramp(phase,empty_region,-1);  
    end
    

end
