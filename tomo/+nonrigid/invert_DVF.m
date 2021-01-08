% INVERT_DVF simple iterative method for estimation of the inverse deformation
% field 
% Chen, Mingli, et al. "A simple fixed‐point approach to invert a deformation field a." Medical physics 35.1 (2008): 81-88.
% 
% [deform_tensors,inv_deform_tensors] = invert_DVF(deform_tensors, Npix_vol)
% 
%    **deform_tensors    cell of 3D arrays containing forward deformation DVF       
%    **Npix_vol          pixel size of the recosntructed volume      
% Outputs: 
%    ++deform_tensors_linear        calculate forward DVF 
%    ++inv_deform_tensors_linear    calculate inverse DVF 

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2018 by Paul Scherrer Institute (http://www.psi.ch)    |
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

function [deform_tensors,inv_deform_tensors] = invert_DVF(deform_tensors, Npix_vol)


   Niter = 10; 
   %% find invert transformation 
   N = length(deform_tensors); 
    for block = 1:N
        % init guess 
        for ax = 1:3
            inv_deform_tensors{block}{ax} = -deform_tensors{block}{ax}; 
        end
        for i = 1:Niter
            for ax = 1:3
                scale = size(deform_tensors{block}{ax}, ax) / Npix_vol(ax) ; % calculate the deformation in deform_tensors grid 
                inv_deform_tensors{block}{ax} =  inv_deform_tensors{block}{ax}*0.5 + 0.5*utils.interp3_gpu(-deform_tensors{block}{ax}, ...
                                    scale*inv_deform_tensors{block}{1},scale*inv_deform_tensors{block}{2},scale*inv_deform_tensors{block}{3});
            end
        end
        % in my code I assume that deform_tensors and inv_deform_tensors
        % have the same direction (in the astra the direction is swapped)
        for ax = 1:3
            inv_deform_tensors{block}{ax} = -gather(inv_deform_tensors{block}{ax}); 
        end
    end
 
end
