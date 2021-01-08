%HDF5_MV_DATA Move data within an HDF5 file
% hdf5_mv_data creates a new (UNIX-like) hard link at loc_dest to the dataset at
% loc_origin and deletes the hard link to the dataset at loc_origin.
%
% file...       HDF filename
% loc_origin... location of the data that needs to be moved
% loc_dest...   destination and name of the new data
%
% EXAMPLE:
%   % move dataset probe from root to group measurements
%   hdf5_mv_data('./awesome_file.h5', 'probe', 'measurements/probe')
%
%   Please notice that all groups and datasets have to exist before running
%   the script!

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
function hdf5_mv_data( file, loc_origin, loc_dest)
import io.HDF.*
plist = 'H5P_DEFAULT';

fileID = H5F.open(file,'H5F_ACC_RDWR',plist);

gpath1 = strsplit(rm_delimiter(loc_origin), '/');
gid1 = add_groups(fileID, gpath1(1:end-1), plist, false);

gpath2 = strsplit(rm_delimiter(loc_dest), '/');
gid2 = add_groups(fileID, gpath2(1:end-1), plist, false);

try
    datasetID = H5D.open(gid1{end}, gpath1{end});
    dataset = true;
catch
    gid1 = add_groups(fileID,gpath1,plist, false);
    datasetID = gid1{end};
    dataset = false;
end

H5O.link(datasetID,gid2{end},gpath2{end},plist,plist);
if dataset
    H5L.delete(gid1{end}, gpath1{end}, plist);
else
    H5L.delete(gid1{end-1}, gpath1{end}, plist);
end

H5F.close(fileID);


end
