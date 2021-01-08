% DECOMPOSE_AFFINE_MATRIX calculate rotation, shear, asymmetry and scale from affine matrix 
% 
%  [scale, asymmetry, rotation, shear]  = decompose_affine_matrix(affine_mat)
% 
% Inputs:
% ** affine_mat  affine matrix = A1*A2*A3*A4
%
% returns: 
% ++scale      A1 = [scale, 0; 0, scale]
% ++asymmetry   A2 = [1+asymmetry/2,0; 0,1-asymmetry/2]
% ++rotation    A3 = [cosd(rotation), sind(rotation); -sind(rotation), cosd(rotation)]
% ++shear       A4 = [1,0;tand(shear),1];
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



function [scale, asymmetry, rotation, shear]  = decompose_affine_matrix(affine_mat)
        err_fun = @(x)(affine_mat - x(1)*[1+x(2)/2,0; 0,1-x(2)/2]*[cosd(x(3)), sind(x(3)); -sind(x(3)), cosd(x(3))] * [1,0;tand(x(4)),1]);
        options = optimoptions('lsqnonlin','Display','off');
        xopt = lsqnonlin( err_fun, [1,0,0,0],[],[],options); 
        scale = xopt(1);
        asymmetry = xopt(2);
        rotation = xopt(3);
        shear = xopt(4); 
end
