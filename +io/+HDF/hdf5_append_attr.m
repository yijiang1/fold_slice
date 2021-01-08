%HDF5_APPEND_ATTR Append an attribute to a dataset
%
% file...       HDF filename
% attr...       structure of attributes
% loc...        location of the dataset that needs to be removed
% 
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

function hdf5_append_attr( file, attr, loc)
import io.HDF.*
plist = 'H5P_DEFAULT';

fileID = H5F.open(file,'H5F_ACC_RDWR',plist);

gpath = strsplit(rm_delimiter(loc), '/');

gpath_depth = length(gpath);
gid{1} = fileID;
for ii=1:gpath_depth
    try
        gid{end+1} = H5G.open(gid{ii}, gpath{ii}, plist);
    catch
        gid{end+1} = H5D.open(gid{ii}, gpath{ii});
    end
end

fn = fieldnames(attr);
for ii=1:length(fn)
    write_attribute(gid{end}, attr.(fn{ii}), fn{ii}, true);
end
    
H5F.close(fileID);


end
