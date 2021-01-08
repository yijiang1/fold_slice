%  UNWRAP2D_FFT2 simple and very fast 2D phase unwrapping based on FT
%  integration of DIC signal (phase gradient) 
%
%  [phase, residues] = unwrap2D_fft2(img, empty_region=[], step=0, weights=[], polyfit_order=1)
% 
% Inputs: 
%   **img           - either complex valued image or real valued phase gradient 
% *optional*
%   **empty_region  - 2x1 or 1x1 vector, size of empty region assumed around edges for phase offset removal 
%   **step          - used to calculate finite difference gradient, 0 = analytical (default) expression
%   **polyfit_order   -1 = dont assume anything about the removed phase,
%                       subtract linear offset from each horizontal line separatelly
%                      0 = assume that removed degree of freedom is only a constant offset
%                      1 = assume that removed degree of freedom is a 2D plane 
% *returns*: 
%  ++phase          - unwrapped phase with subtracted phase ramp / offset 
%  ++residues       - calculated binary map of phase residuas 
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




function [phase, residues] = unwrap2D_fft2(img, empty_region, step, weights, polyfit_order)
    import utils.* 
    import math.* 
    
    if isreal(img)
        error('Complex-valued input array was expected')
    end
    
    if nargin < 3 || isempty(step)
        step = 0; % use analytical method 
    end
    if nargin < 4 || isempty(weights)
        weights = 1; 
    end
    if nargin < 5
        polyfit_order = 1; % subtract 2D plane 
    end
    
    weights = max(0, min(1, single(weights)));  % weights are assumed in range [0,1]
    
    img = weights .* img ./ (abs(img)+eps);
    clear weights 
    
    padding = [64,64];  % padding to avoid periodic boundary artefacts 
    if any(padding > 0)
        img = padarray(img,padding,'symmetric','both'); 
        img = smooth_edges(img, 5, [1,2]);
    end

    
    [dX,dY] = get_phase_gradient_2D(img, step, 0);
    clear img 
        
    phase = real(get_img_int_2D(dX,dY));
    
    if any(padding > 0)
        % remove padding 
        ind = {padding(1):size(phase,1)-padding(1)-1,padding(2):size(phase,2)-padding(2)-1, ':'};
        phase = phase(ind{:}); 
    end
    
    if exist('empty_region', 'var') && ~isempty(empty_region)
       phase = remove_sinogram_ramp(phase,empty_region, polyfit_order);  
    end
    
    if nargout > 1
        residues = abs(findresidues(img)) > 0.1; 
    end

end
