% SVD_regularize_3D_fields - use SVD to constraint the reconstructed DVF
% and enforce smoothness in the DVF reconstruction 
%
%   shift_3D_total = SVD_regularize_3D_fields(shift_3D_total, Nsvd, SVD_smoothing)
%
% Inputs:
%    **shift_3D_total - (cell of 3D arrays) recovered deformation field 
%    **Nsvd - (int), maximal number of SVD modes to be recovered (should be less than number of subtomos)
%    **SVD_smoothing - (scalar), smoothing constant between 0 to 0.25
% returns:
%    ++shift_3D_total - (cell of 3D arrays) regularized deformation field 


%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
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


function shift_3D_total = SVD_regularize_3D_fields(shift_3D_total, Nsvd, SVD_smoothing)

        % SVD regularization 
        Nblocks= length(shift_3D_total);
        Nps = size(shift_3D_total{1}{1});
        for kk = 1:3
            for ll = 1:Nblocks
                shift_3D_mat(:,:,:,kk,ll) = shift_3D_total{ll}{kk};
            end
        end
        shift_3D_mat = reshape(shift_3D_mat,[],Nblocks);
        [U,S,V] = fsvd(shift_3D_mat, min(Nsvd, Nblocks));
        %% apply a bit of smoothness 
        kernel = [SVD_smoothing, 1-2*SVD_smoothing, SVD_smoothing]'; 
        V = conv2(V, kernel, 'same') ./ conv2(ones(Nblocks,min(Nsvd, Nblocks)), kernel, 'same') ; 
        shift_3D_mat = U*S*V';
        shift_3D_mat = reshape(shift_3D_mat, [Nps,3,Nblocks]);
        for kk = 1:3
            for ll = 1:Nblocks
               shift_3D_total{ll}{kk} = shift_3D_mat(:,:,:,kk,ll);
            end
        end
        
end