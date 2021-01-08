% MAKE_SYNTHETIC_PROJECTIONS   creates projections that can be used as
% initial guess for ptychography reconstruction in order to solver
% iterativelly the ptychotomo task 
%
% merged = make_synthetic_projections(stack_object, sinogram_abs, sinogram_phase,total_shift,object_ROI)
%
%   Inputs: 
%       **stack_object - measured projections 
%       **sinogram_abs - aligned absorbtion sinogram 
%       **sinogram_phase  - aligned phase sinogram 
%       **total_shift - reconstructed shifts oft the measured projections 
%       **object_ROI - reliable region of the projections 
%   Outputs: 
%       merged - merged projections using stack_object, sinogram_abs and sinogram_phase


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
%   in written form in the publication: "Data processing was carried out 
%   using the "cSAXS matlab package" developed by the CXS group,
%   Paul Scherrer Institut, Switzerland." 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided "as they are" without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.


function merged = make_synthetic_projections(stack_object, sinogram_abs, sinogram_phase,total_shift,object_ROI)

    [Nlayers,Nw,~] = size(sinogram_abs); 
    Npx_proj = [size(stack_object,1),size(stack_object,2)]; 

    assert(length(object_ROI{1}) == Nlayers &&  length(object_ROI{2}) == Nw, 'Reconstructed tomograms have to contain the full field of view, ie vert_range = object_ROI{1}' )
    
    if isreal(sinogram_abs) && isreal(sinogram_phase)
        sinogram_abs = exp(-sinogram_abs); 
        sinogram_phase = exp(-1i*sinogram_phase); 
    else
        sinogram_phase = sinogram_abs ./ (abs(sinogram_abs)+1e-3); 
        sinogram_abs = abs(sinogram_abs); 
    end
    

    % find reliability region 
    win = tukeywin(Nw, 0.2)'.*tukeywin(Nlayers, 0.2); 
    win = utils.crop_pad(win, Npx_proj);
    % find amplitude correction 
    corr = math.sum2(abs(stack_object(object_ROI{:},:)) .* sinogram_abs) ./ math.sum2(sinogram_abs.^2); 
    % merge phase and amplitude from tomo and measurements 
    merged =  win.*utils.crop_pad(sinogram_phase,Npx_proj) + (1-win).*stack_object ./ (abs(stack_object)+1e-2); 
    % enforce phasor 
    merged = merged ./ (abs(merged) + 1e-2); 
    % apply amplitude 
    merged = merged .* ((win .* corr .* utils.crop_pad(sinogram_abs,Npx_proj) + (1-win).* abs(stack_object)));
    
    clear sinogram_phase sinogram_abs
    
    if ~isempty(total_shift)
        % apply shift to match the original data 
        merged = utils.imshift_fft(merged, -total_shift); 
    end

end
