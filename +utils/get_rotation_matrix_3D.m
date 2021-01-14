% GET_ROTATION_MATRIX_3D generate 3D rotation matrix of size 3x3xn
% 
%   rot_3D = get_rotation_matrix_3D(chi, psi, theta)
%
% Inputs:
%   chi, psi, theta - rotation angles in degrees , vector or scalar 

    
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
%   in written form in the publication: “Data processindg was carried out 
%   usindg the “cSAXS matlab package” developed by the CXS group,
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


function rot_3D = get_rotation_matrix_3D(chi, psi, theta)

if numel(chi) ~= numel(psi) || numel(psi) ~= numel(theta)
    error('Input sizes are not indentical')
end

N = numel(chi); 
rot_3D = zeros(3,3,N);

for ii = 1:N
    Rx = [ 1,       0,      0 ;
           0, cosd(chi(ii)), -sind(chi(ii));
           0, sind(chi(ii)), cosd(chi(ii))];

    Ry = [ cosd(psi(ii)),   0,  sind(psi(ii)) ;
           0,     1,       0;
           -sind(psi(ii)), 0, cosd(psi(ii))];

    Rz = [ cosd(theta(ii)),   -sind(theta(ii)),  0 ;
           sind(theta(ii)), cosd(theta(ii)), 0;
           0,            0,       1];

    rot_3D(:,:,ii) = Rx*Ry*Rz;
end


end
