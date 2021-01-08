%WRITE_ATTRIBUTE write attribute data_name with value data to ID gid

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

function write_attribute(gid, data, data_name, varargin)
import io.HDF.*

if nargin > 3
    safe = varargin{1};
else
    safe = false;
end

[datatypeID, data] = get_datatype(data);
if ischar(data)
    % The ptycho C++ code expects strings as H5S_SCALAR, so we have to
    % convert it
    data = {data};
    filetype = H5T.copy ('H5T_FORTRAN_S1');
    H5T.set_size (filetype,'H5T_VARIABLE');
    memtype = H5T.copy ('H5T_C_S1');
    H5T.set_size (memtype, 'H5T_VARIABLE');
    space = H5S.create ('H5S_SCALAR');
    if safe
        try
            attr = H5A.create (gid, data_name, filetype, space, 'H5P_DEFAULT');
        catch
            H5A.delete(gid, data_name);
            attr = H5A.create (gid, data_name, filetype, space, 'H5P_DEFAULT');
        end
    else
        attr = H5A.create (gid, data_name, filetype, space, 'H5P_DEFAULT');
    end

    H5A.write (attr, memtype, data);
    
elseif iscell(data)
    % If it is a cell, save it as 1D dataset 
    H5T.set_size(datatypeID,'H5T_VARIABLE');    
    agcv = H5ML.get_constant_value('H5S_UNLIMITED');
    dspace = H5S.create_simple(1,numel(data),agcv);
    
    plist = H5P.create('H5P_ATTRIBUTE_CREATE');
    if safe
        try
            attr = H5A.create(gid,data_name,datatypeID,dspace,plist);
        catch
            H5A.delete(gid, data_name);
            attr = H5A.create(gid,data_name,datatypeID,dspace,plist);
        end
    else
        attr = H5A.create(gid,data_name,datatypeID,dspace,plist);
    end
    H5A.write(attr,'H5ML_DEFAULT',data);
    
else 
    
    acpl = H5P.create('H5P_ATTRIBUTE_CREATE');
    dims = size(data);    
    if length(dims)>1 && dims(2)~=1
        if dims(1) == 1
            space_id = H5S.create_simple(dims(1), dims(2), []);
        else
            space_id = H5S.create_simple(dims(1), dims, []);
        end
    else
        space_id = H5S.create('H5S_SCALAR');
    end    
    if safe
        try
            attr = H5A.create(gid,data_name,datatypeID,space_id,acpl);
        catch
            H5A.delete(gid, data_name);
            attr = H5A.create(gid,data_name,datatypeID,space_id,acpl);
        end
    else
        attr = H5A.create(gid,data_name,datatypeID,space_id,acpl);
    end


    H5A.write(attr,'H5ML_DEFAULT',data)
    

end
H5A.close(attr);
end

