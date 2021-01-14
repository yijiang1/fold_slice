% MTIMES_STACK  Extension of @mtimes function for stacked images and GPU
%
%    C = mtimes_stack(A,B)
%     returns the propagated wavefield
%    Inputs: 
%       **A           first matrix
%       **B           second matrix 
%    *returns*
%       ++C           product matrix 
%
    
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


function C = mtimes_stack(A,B)

    if isa(A, 'gpuArray') || isa(B, 'gpuArray')
        C = pagefun(@mtimes, A,B); 
        return
    end

    % CPU code 
    if ismatrix(A) && ismatrix(B)
        C = mtimes(A,B); 
    elseif ismatrix(A) && ndims(B) == 3
        Np = size(B);    
        C = reshape(A*reshape(B,Np(1),[]), Np);
    elseif ndims(A) == 3 && ismatrix(B)
        Np = size(A);    
        fdims = 1:ndims(A);
        fdims(1:2) = [2,1]; 
        C = permute(reshape(B*reshape(permute(A,fdims),Np(2),[]),Np(fdims)),fdims);
    else
        C = zeros(size(A,1), size(B,2), size(B,3), 'like', A);
        for ii = 1:size(B,3)
           C(:,:,ii) = A(:,:,ii) * B(:,:,ii); 
        end
    end

end

