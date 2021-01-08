% Determine datatype for HDF files

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

function [datatype_h5, data] = get_datatype(data)

if isempty(data)
    data = '';
end

switch class(data)
    case 'uint32'
        datatype_h5 = 'H5T_STD_U32LE';
    case 'int32'
        datatype_h5 = 'H5T_STD_I32LE';
    case 'int64'
        datatype_h5 = 'H5T_STD_I64LE';
    case 'uint64'
        datatype_h5 = 'H5T_STD_U64LE';
    case 'double'
        if isreal(data)
            datatype_h5 = 'H5T_NATIVE_DOUBLE';
        else
            datatype_h5 = 'complex';
        end
    case 'single'
        if isreal(data)
            datatype_h5 = 'H5T_NATIVE_FLOAT';
        else
            datatype_h5 = 'complex';
        end
    case 'logical'
        data = uint32(data);
        datatype_h5 = 'H5T_STD_U32LE';
    case 'char'
        if size(data,1)>1
            data = cellstr(data);
            datatype_h5 = 'char_array';
        else
            datatype_h5 = 'H5T_UNIX_D32BE'; %'H5T_UNIX_D32LE';
        end
    case 'cell'
        datatype_h5 = 'H5T_C_S1';
    case 'struct'
        datatype_h5 = 'struct';
        
    otherwise
        error('Unknown data type %s', class(data))
end

end

