% Read interferometer positions written by Orchestra
% Input is the filename with path
% Output is a structure containing fields:
%    The two values of the one line header originally 'Scan' and 'Samroy'
%    Values for each point of 10 expected columns of numbers
% 12 June 2013
% June6 2015 - Changed in order to accept an arbitrary number
% of values in order to be compatible with 10 columns for flOMNI and 19 for
% OMNY
% This function should be replaced by beamline.read_position_file in the 
% ptycho codes and deprecated.

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

function struct_out = read_omny_pos( omnyposfile )
%disp(omnyposfile)
struct_out = beamline.read_position_file( omnyposfile );
%disp(size(struct_out.TotalPoints))
end

