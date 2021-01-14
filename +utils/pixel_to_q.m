% Call function without arguments for instructions on how to use it

% Filename: $RCSfile: pixel_to_q.m,v $
%
% $Revision: 1.1 $  $Date: 2008/06/10 17:05:14 $
% $Author:  $
% $Tag: $
%
% Description:
% calculated momentum transfer q in inverse Angstroem from pixel numbers
% relative to the beam center
%
% Dependencies: 
% none
%
% history:
%
% June 9th 2008: 1st documented version

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

function [ q_A ] = pixel_to_q( pixel, pixel_size_mm, det_dist_mm, E_keV )

if (nargin ~= 4)
    fprintf('Usage:\n');
    fprintf('[ q_A ] = %s( pixel, pixel_size_mm, det_dist_mm, E_keV );\n',...
        mfilename);
    error('Wrong number of parameters, 4 expected, %d found',nargin);
end

lambda_A = 12.39852 / E_keV;

q_A = 4*pi * sin( atan(pixel*pixel_size_mm/det_dist_mm) /2) / lambda_A;
