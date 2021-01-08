%WRITE_DATASET write dataset data_name, containing data to ID gid

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

function write_dataset(data, gid, data_name, plist, varargin)
import io.HDF.*

extend_data = false;
link = false;
link_type = '';
cellstrdata = false;
write_data = true;

if ~isempty(varargin)
    comp = varargin{1};
else
    comp = true;
end
if nargin > 5
    overwrite = varargin{2};
else
    overwrite = true;
end

if nargin > 6
    data_attr = varargin{3};
else
    data_attr = [];
end

if nargin > 7
    extend_dim = varargin{4};
else
    extend_dim = 0;
end

if nargin > 8
    extendable = varargin{5};
else
    extendable = false;
end

if nargin > 9
    extend_offset = varargin{6};
else
    extend_offset = 0;
end

if nargin > 10
    extend_maxdims = varargin{7};
else
    extend_maxdims = 0;
end

[datatype, data] = get_datatype(data);
filespaceID = [];

    function create_dataspace()
        
        if extendable
            unlimited = H5ML.get_constant_value('H5S_UNLIMITED');
            dims_max = repmat(unlimited, 1, numel(dims));
        else
            dims_max = dims;
        end
        
        if ~extend_dim
            if ~extendable
                dataspaceID = H5S.create_simple(length(dims), fliplr(dims), fliplr(dims_max));
            else
                try 
                    datasetID = H5D.open(gid, data_name);
                    filespaceID = H5D.get_space(datasetID);
                    [~, spaceDims] = H5S.get_simple_extent_dims(filespaceID);
%                     spaceDims = fliplr(spaceDims);
                    
                    start = ones(1,numel(dims))-1;
                    count = dims;

                    stride = ones(1, numel(start));
                    boundsEnd = start + (count).*stride;
                    new_dims = fliplr(boundsEnd);
                    H5S.close(filespaceID);
                    H5D.set_extent(datasetID,new_dims);
                    filespaceID = H5D.get_space(datasetID);
                    H5S.select_hyperslab(filespaceID, 'H5S_SELECT_SET', fliplr(start), fliplr(stride), ...
                        fliplr(count), ones(1,length(start)));
                    
                    dataspaceID = H5S.create_simple(numel(count),fliplr(count),[]);
                    
                    extend_data = true;

                catch
                    dataspaceID = H5S.create_simple(length(dims), fliplr(dims), fliplr(dims_max));
                end
            end

                    
        else
            try
                datasetID = H5D.open(gid, data_name);
                filespaceID = H5D.get_space(datasetID);
                [~, spaceDims] = H5S.get_simple_extent_dims(filespaceID);
                spaceDims = fliplr(spaceDims);
                if extend_offset
                    start = [ones(1,extend_dim-1) extend_offset+1]-1;
                else
                    start = [ones(1,extend_dim-1) spaceDims(end)+1]-1;
                end

                if numel(spaceDims) > numel(dims)
                    count = [dims 1];
                else
                    count = dims;
                end
                stride = ones(1, numel(start));
                boundsEnd = start + (count-1).*stride;
                if extend_maxdims
                    boundsStart = spaceDims;
                    boundsStart(end) = extend_maxdims;
                else
                    boundsStart = spaceDims;
                end
                new_dims = fliplr(max(boundsStart,boundsEnd+1));
                H5S.close(filespaceID);
                H5D.set_extent(datasetID,new_dims);
                filespaceID = H5D.get_space(datasetID);
                H5S.select_hyperslab(filespaceID, 'H5S_SELECT_SET', fliplr(start), fliplr(stride), ...
                    fliplr(count), ones(1,length(start)));
                
                dataspaceID = H5S.create_simple(numel(count),fliplr(count),[]);
                
                extend_data = true;
                              
                
            catch
                unlimited = H5ML.get_constant_value('H5S_UNLIMITED');
                maxdims = [dims(1:extend_dim-1) unlimited];
%                 maxdims = repmat(-1, 1, extend_dim);
                if numel(maxdims) > numel(dims)
                    dims = [dims 1];
                end
                dataspaceID = H5S.create_simple(length(dims), [fliplr(dims)], fliplr(maxdims));
            end
        end
    end

if strcmp(datatype, 'complex')
    
    %%% prepare compound dataset for complex input data
    dims = size(data);

    data_temp = data;
    data = [];
    data.r = real(data_temp);
    data.i = imag(data_temp);
    
    create_dataspace();
    
    % Create the required data types
    complexType = H5T.copy(get_datatype(data.r));
    sz = H5T.get_size(complexType);
    
    % Create the compound datatype for memory.   
    datatypeID = H5T.create ('H5T_COMPOUND', 2*sz);
    H5T.insert (datatypeID, 'r',0, complexType);
    H5T.insert (datatypeID, 'i',sz, complexType);
    memtype = datatypeID;
    
    data_attr.MATLAB_class =  'complex';

    
elseif ischar(data)
    % check if char is a link
    ch_entrs = strsplit(data, ':');
    if length(ch_entrs) >= 2
        link = true;
        if strcmp(ch_entrs{1}, 'ext')
            % prepare external link
            link_type = 'ext';
        elseif strcmp(ch_entrs{1}, 'int_soft')
            % prepare internal soft link
            link_type = 'int_soft';
        elseif strcmp(ch_entrs{1}, 'int_hard')
            % prepare internal hard link
            link_type = 'int_hard';
        end
    else
        data = {data};
        datatypeID = H5T.copy ('H5T_FORTRAN_S1');
        H5T.set_size (datatypeID,'H5T_VARIABLE');
        memtype = H5T.copy ('H5T_C_S1');
        H5T.set_size (memtype, 'H5T_VARIABLE');
        dataspaceID = H5S.create ('H5S_SCALAR');

    end
    
elseif iscell(data) || strcmp(datatype, 'char_array')

    
    if iscellstr(data)
        cellstrdata = true;
        datatypeID = H5T.copy ('H5T_C_S1');
        H5T.set_size (datatypeID, 'H5T_VARIABLE');
        
        dgcv = H5ML.get_constant_value('H5S_UNLIMITED');
        dataspaceID    = H5S.create_simple(1,numel(data),dgcv);
        memtype = datatypeID;
        plist_cr = H5P.create('H5P_DATASET_CREATE');
        H5P.set_chunk(plist_cr,1);
        if strcmp(datatype, 'char_array')
            data_attr.MATLAB_class =  'char_array';
        end
    else
        write_data = false;
        fn_names = cell(1,length(data));
        cell_gid = add_groups(gid, data_name, plist, true);
        for ii=1:length(data)
            fn_names{ii} = sprintf([data_name '_%d'],ii-1);
            write_dataset(data{ii}, cell_gid, fn_names{ii}, plist, comp, overwrite);
        end
        write_attribute(cell_gid, 'cell', 'MATLAB_class');


    end
    
elseif isstruct(data)
    
    write_data = false;
    struct_gid = add_groups(gid, data_name, plist, true);
    add_content(data, struct_gid, plist, comp, overwrite);
    
else
    datatypeID = H5T.copy(datatype);
    dims = size(data);
    if isfield(data_attr, 'save2hdf5DataShape')
        dims = data_attr.save2hdf5DataShape;
    end

    % prepare dataspace
    create_dataspace();

    memtype = 'H5ML_DEFAULT';
end

%%% create groups and write data
if comp && ~iscell(data) && ~ischar(data) && write_data || extend_dim || extendable
    % define compression and chunk size
    plist_ch = H5P.create('H5P_DATASET_CREATE');
    if length(dims)>=3
        chunk_dims = [dims(1) dims(2) ones(1, numel(dims)-2)];
    else
        chunk_dims = dims;
    end

    h5_chunk_dims = fliplr(chunk_dims);
    H5P.set_chunk(plist_ch,h5_chunk_dims);
    H5P.set_shuffle(plist_ch);
    if comp
        H5P.set_deflate(plist_ch,comp);
    end
    
    % Try to create a new dataset. If it exists, try to open it.
    try
        if ~extend_data
            if cellstrdata
                datasetID = H5D.create(gid,data_name,datatypeID,dataspaceID,plist_cr);
            else
                datasetID = H5D.create(gid,data_name,datatypeID,dataspaceID,plist_ch);
%                 create_dataspace();
            end
        end
    catch
        if ~overwrite
            try
                datasetID = H5D.open(gid, data_name);
            catch
                error('Could not create dataset %s! Try a different name or overwrite the already existing file.', data_name);
            end
        else
            keyboard
            error('Dataset %s already exists! Try a different name or overwrite the already existing file.', data_name);
        end
    end
    
elseif ~link && write_data
    % Same as above but without compression:
    % Try to create a new dataset. If it exists, try to open it.
    try
        if cellstrdata
            datasetID = H5D.create(gid,data_name,datatypeID,dataspaceID,plist_cr);
        else
            datasetID = H5D.create(gid,data_name,datatypeID,dataspaceID,plist);
        end
    catch
        if ~overwrite
            try
                datasetID = H5D.open(gid, data_name);
            catch
                error('Could not open dataset %s! Try a different name or overwrite the already existing file.', data_name);
            end
        else
            error('Dataset %s already exists! Try a different name or overwrite the already existing file.', data_name);
        end
    end
end

if write_data
    % write data to disk or link data
    
    if ~link && ~extend_data
        H5D.write(datasetID,memtype,'H5S_ALL','H5S_ALL',plist ,data);
        % append attributes if needed
        if ~isempty(data_attr)
            fn = fieldnames(data_attr);
            for ii=1:length(fn)
                write_attribute(datasetID, data_attr.(fn{ii}), fn{ii}, true);
            end
        end
        H5D.close(datasetID);
    elseif extend_data 
        H5D.write(datasetID,memtype,dataspaceID, filespaceID, plist, data)
    elseif strcmp(link_type, 'ext')
        H5L.create_external(ch_entrs{2},ch_entrs{3},gid,data_name,plist,plist);
    elseif strcmp(link_type, 'int_hard')
        error('Currently not supported, sorry!')
        %     H5L.create_hard(ch_entrs{2},'g3',gid1,'g4',plist,plist);
    elseif strcmp(link_type, 'int_soft')
        H5L.create_soft(ch_entrs{2},gid,data_name,plist,plist);
    end
    
    


end

end

