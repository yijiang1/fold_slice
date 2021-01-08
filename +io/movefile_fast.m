% MOVEFILE_FAST  faster alternative to the matlab movefile function 
% 
% movefile_fast(source, destination)
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

function movefile_fast(source, destination)

     Sd = java.io.File(destination); 
     if Sd.isDirectory
         % if provided path to move is a directory, create full file path from the source name 
         [~, filename, ext] = fileparts(source); 
         destination = fullfile(destination, [filename, ext]);
         Sd = java.io.File(destination); 
     end
     Ss = java.io.File(source);
     assert(Ss.canRead, sprintf('File %s does not exist or is not readable', source))
     % move the file using java
     Ss.renameTo(Sd);
     assert(Sd.canWrite, sprintf('Moving file %s to %s failed', source, destination))
end


