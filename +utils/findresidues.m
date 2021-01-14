% function residues = findresidues(phase)
% Receives phase in radians, returns map of residues
% Manuel Guizar - Sept 27, 2011
% R. M. Goldstein, H. A. Zebker and C. L. Werner, Radio Science 23, 713-720
% (1988).
% Inputs
% phase      Phase in radians
% disp      = 0, No feedback
%           = 1, Text feedback (additional computation)
%           = 2, Text and graphic display (additional computation)
% Outputs
% residues  Map of residues, note they are valued +1 or -1

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


function residues = findresidues(phase)

if ~isreal(phase) 
    phase = angle(phase); 
end

residues =            wrapToPi(phase(2:end,1:end-1,:)   - phase(1:end-1,1:end-1,:));
residues = residues + wrapToPi(phase(2:end,2:end,:)     - phase(2:end,1:end-1,:));
residues = residues + wrapToPi(phase(1:end-1,2:end,:)   - phase(2:end,2:end,:));
residues = residues + wrapToPi(phase(1:end-1,1:end-1,:) - phase(1:end-1,2:end,:));
residues = residues/(2*pi);

end

function x = wrapToPi(x)
    x = mod(x+pi, 2*pi)-pi;
end
