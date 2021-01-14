% READ_POSITION_FILE Read positions from a file, the format and header are 
% compatible with  multiple interferometer positions and standard deviations  
% as written by Orchestra and the sgalil spec macro.
%
% struct_out = read_position_file( posfile )
% Inputs:
%   **posfile   filename with path
% *returns*
%   ++struct_out    is a structure containing fields including the values
%   for header and for each scanning point

% 12 June 2013
% June6 2015 - Changed in order to accept an arbitrary number
% of values in order to be compatible with 10 columns for flOMNI and 19 for
% OMNY
% 15 Apr 2019 - Changed name and generalized description beyond OMNY and
% Orchestra

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

function struct_out = read_position_file( posfile )

assert(exist(posfile, 'file')>0, ['Position file ', posfile, ' not found'])

f = fopen(posfile,'r');
header = textscan(f,'%s %d, %s %f',1);
struct_out.(header{1}{1}) = header{2};
struct_out.(header{3}{1}) = header{4};
names =  textscan(f,'%s',1,'Delimiter','\r');
names = strsplit(char(names{1}));
reading_string =  ['%f', repmat(' %f',1,numel(names)-1)];
values = textscan(f,reading_string);
fclose(f);

for ii = 1:numel(names)
    struct_out.(char(names(ii))) = values{ii};
end

end

