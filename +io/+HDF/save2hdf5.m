%SAVE2HDF5 saves matlab data to a Hierarchical Data Format file (hdf5)
%   
%   filename...     full path to file, including file extension
%   data...         matlab structure or array or link
%   data_name...    needed if input data is not a matlab structure, needs
%                   to be given as name/value pair
%
%   *optional*
%   overwrite...    replace existing file if it exists
%   gpath...        specify the group to which you want to append the data
%                   (only if data is an array); default root ('/')
%   Attributes...   structure of attributes; will be appended to current
%                   gpath
%   comp...         compression level; default 0 (no compression)
%   creator...      attribute in root; default 'ptycho_recons'
%   
%
%   If you want to save a structure, everything declared within an 'Attributes'
%   fieldname will be treated as an attribute to the current group.
%   If you want to add attributes to a dataset, you have to define your
%   data within .Value and your attributes within .Attributes.
%
%   A simple structure could look like:
%       h5_struc = [];
%       h5_struc.probe_mask = ones(256,256);
%       h5_struc.Attributes.probe_id = 1;
%       h5_struc.measurement.n0.diff = fmag(:,:,1);
%       h5_struc.measurement.n0.Attributes.detector = 0;
%       h5_struc.measurement.n1.diff.Value = fmag(:,:,2);
%       h5_struc.measurement.n1.diff.Attributes.slice = 2;
%
%   fmag(:,:,1) will be written to dataset 'diff' in group '/measurement/n0' 
%   fmag(:,:,2) with attribute 'slice' will be written to dataset 'diff' in 
%   group '/measurement/n1' 
%
%
%   EXAMPLES:
%     -) if data is a matlab structure:
%           save2hdf5('./awesome_file.h5', data);
%           save2hdf5('./awesome_file.h5', data, 'overwrite', true);
%                    
%
%     -) if data is a matlab array:
%           save2hdf5('./awesome_file.h5', data, 'data_name', data_name);
%           save2hdf5('./awesome_file.h5', data, 'data_name', 'my_dataset',...
%                     'gpath', 'group1/group2', 'Attributes', attr_struc);
%
%     -) if data is a link:
%           currently, only external links ('ext') and internal soft links
%           ('int_soft') are supported
%
%           external links have to be specified by a single string with 
%           3 sections: '<link_type>:<file_path>:<target_object>'
%               
%               e.g.: 'ext:./awesome_file2.h5:/data' 
%                   save2hdf5('./awesome_file.h5',...
%                     'ext:./awesome_file2.h5:/data', 'data_name', data_name)
%               
%               will create a link called $data_name to dataset (or group) '/data' 
%               in './awesome_file2.h5'
%
%           internal links have to be specified by a single string with
%           2 sections: '<link_type>:<target_object>'
%
%               e.g.: 'int_soft:/data' 
%                   save2hdf5('./awesome_file.h5',...
%                     'int_soft:/data', 'data_name', data_name, 'gpath', 'g1/g2')
%               
%               will create a link called $data_name to dataset (or group) '/data' 
%               in '/g1/g2'
%
%
%    Please notice that structures are not supported as attributes, i.e.
%       h5_struc = [];
%       h5_struc.attr.probe.probe_id = 1;
%
%       save2hdf5('./awesome_file.h5', h5_struc)
%
%    will crash! 
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

function save2hdf5( filename, data, varargin)
import io.HDF.*

% take care of input arguments
overwrite = false;
gpath_full = '';
attr = [];
data_name = '';
comp = 0;
creator = 'ptycho_recons';
iscopy = false;
extend_dim = 0;
extendable = false;
extend_offset = 0;
extend_maxdims = 0;

vararg = cell(0,0);
% parse the variable input arguments vararg = cell(0,0);
if ~isempty(varargin)
    for ind = 1:2:length(varargin)
        name = varargin{ind};
        value = varargin{ind+1};
        switch lower(name)
            case 'data_name'
                data_name = value;
            case 'overwrite'
                overwrite = value;
            case 'gpath'
                gpath_full = value;
            case 'attr'
                attr = value;
            case 'comp'
                comp = value;
            case 'creator'
                creator = value;
            case 'iscopy'
                iscopy = value;
            case 'extend_dim'
                extend_dim = value;
            case 'extendable'
                extendable = value;
            case 'extend_offset'
                extend_offset = value;
            case 'extend_maxdims'
                extend_maxdims = value;
                
            otherwise
                vararg{end+1} = name;
                vararg{end+1} = value;
        end
    end
end

if ~isstruct(data)
    full_data = false;
else
    full_data = true;
end

if ~isstruct(data) && isempty(data_name)
    data_name = inputname(2);
    if isempty(data_name)
        error('Please specify the data_name.')
    end
end

if extendable && extend_dim
    error('Extending the dimension of an unlimited dataset is currently not supported.');
end

plist = 'H5P_DEFAULT';

%%% create file if it does not exist
if exist(filename, 'file')&&~overwrite
    fileID = H5F.open(filename,'H5F_ACC_RDWR',plist);
else
    fileID = H5F.create(filename,'H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT');
    if ~iscopy
        write_attribute(fileID, filename, 'filename');
        write_attribute(fileID, creator,'creator');
        write_attribute(fileID, datestr(now),'file_time');
    end
end

if full_data
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%% data as structure %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    add_content(data, fileID, plist, comp, overwrite)
       
    
else
    
    %%%%%%%%%%%%%%%%%%%%%
    %%% data as array %%%
    %%%%%%%%%%%%%%%%%%%%%

    % prepare group handles
    if ~isempty(gpath_full)
        gpath = strsplit(rm_delimiter(gpath_full), '/');
        gid = add_groups(fileID, gpath, plist, false);
    else
        gid{1} = fileID;
    end
    
    % write data to file
    write_dataset(data, gid{end}, data_name, plist, comp, overwrite, [], extend_dim, extendable, extend_offset, extend_maxdims);
    
    % append attributes
    if ~isempty(attr)
        attr_fn = fieldnames(attr);
        for ii=1:length(attr_fn)
            write_attribute(gid{end}, attr.(attr_fn{ii}), attr_fn{ii}, true);
        end
    end
    
end

% close handles
H5F.close(fileID);

end


