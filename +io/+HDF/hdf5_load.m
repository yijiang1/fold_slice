% HDF5_LOAD Load an hdf5 file
%
% DATA = HDF5_LOAD(filename) reads a complete file hierarchy recursively, with
% file name/path being specified by the 'filename' argument
%
% DATA = HDF5_LOAD(filename, '-a') reads a complete file hierarchy
% recursively, including attributes
%
% DATA = HDF5_LOAD(filename, location) reads a particular group, link, or a single dataset
% specified by the 'location' argument
%
% ATT = HDF5_LOAD(filename, location, '-a') reads all datasets and attributes associated
% with a particular location in the file (group, link or dataset)
%
% ATT = HDF5_LOAD(filename, location, '-ca') reads all datasets and attributes associated
% with a particular location in the file (group, link or dataset) and
% converts datasets to a specific matlab class based on attribute 'MATLAB_class'
%
% SLICE = HDF5_LOAD(filename, location, {rowRange, colRange, frameRange, ...}) reads a
% portion of a dataset along specified dimentions, where slicing ranges can be defined in
% the following ways (negative indexes count from the end of the corresponding dimensions):
%   range = scalar_index - reads a particular row/col/frame/... (indentical to
%       'range = [scalar_index, scalar_index]')
%   range = [start_index, end_index] - reads all data between start and end
%       indexes
%   range = [start_index, Inf] - reads all data from start_index to the last
%       existing element in the file
%   range = [], or range is omitted at the end - reads the full range of values for that
%       dimention (indentical to 'range = [1, Inf]')
%
% Examples:
%   hdf5_load('scan_003.hdf5')
%   hdf5_load('scan_003.hdf5', '/entry/sample/description')
%   hdf5_load('scan_003.hdf5', '/entry/collection/data/spec', '-a')
%   hdf5_load('scan_003.hdf5', '/entry/instrument/Pilatus_2M/data', {5})
%   hdf5_load('scan_003.hdf5', '/entry/instrument/Pilatus_2M/data', {[-100, Inf]})
%   hdf5_load('scan_003.hdf5', '/entry/instrument/Pilatus_2M/data', {5, [500, Inf], [1, 100]})
%   hdf5_load('scan_003.hdf5', '/entry/instrument/Pilatus_2M/data', {[], [], [1, 100]})

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
%   using the “cSAXS matlab package” developed by the CXS group
%   and the Science IT group, Paul Scherrer Institut, Switzerland.”
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

function data = hdf5_load(filename, varargin)
import io.HDF.*

load_attr = false;
convert2matlab = false;

narginchk(1, 3);
if nargin == 1
    % Read the complete file hierarchy recursively
    try
        info = h5info(filename);
        info.Name = ''; % a special case of the root group
        
    catch ME
        if strcmp(ME.identifier, 'MATLAB:imagesci:h5info:fileOpenErr')
            ME = MException('hdf5_load:h5info', ...
                strjoin({'File', filename, 'does not exist'}));
        end
        
        throwAsCaller(ME);
    end
    
    [data, links] = hdf5_loadGroup(filename, info);
    data = assign_links(data, info, links);
    
elseif nargin == 2
    if any(strcmp(varargin{1}, {'-a', '-ca', '-c'}))
        % second argument is an attribute flag
        try
            info = h5info(filename);
            info.Name = ''; % a special case of the root group
            
        catch ME
            if strcmp(ME.identifier, 'MATLAB:imagesci:h5info:fileOpenErr')
                ME = MException('hdf5_load:h5info', ...
                    strjoin({'File', filename, 'does not exist or is not valid h5 file'}));
            end
            
            throwAsCaller(ME);
        end
        
        if any(strcmp(varargin{1}, {'-a', '-ca'}))
            load_attr = true;
        end
        
        if any(strcmp(varargin{1}, {'-ca', '-c'}))
            convert2matlab = true;
        end
        
        [data, links] = hdf5_loadGroup(filename, info, convert2matlab, load_attr);
        data = assign_links(data, info, links, varargin{1});
        
    else
        % Read a group or a single dataset
        location = varargin{1};
        try
            info = h5info(filename, location);
            if strcmp(info.Name, '/') % a special case of the root group
                info.Name = '';
            end
            
        catch ME
            if strcmp(ME.identifier, 'MATLAB:imagesci:h5info:fileOpenErr')
                ME = MException('hdf5_load:h5info', ...
                    strjoin({'File', filename, 'does not exist'}));
                
            elseif strcmp(ME.identifier, 'MATLAB:imagesci:h5info:libraryError')
                ME = MException('hdf5_load:h5info', ...
                    strjoin({'H5Location', location, 'was not found in', filename, 'file'}));
            end
            
            throwAsCaller(ME);
        end
        
        if isfield(info, 'Groups')
            % Read a group with its internal hierarchy
            [data, links] = hdf5_loadGroup(filename, info);
            data = assign_links(data, info, links);
            
        elseif isfield(info, 'Datatype')
            % Read a data set
            type = info.Datatype.Class;
            data = hdf5_loadDataset(filename, location, type);
            
        elseif isfield(info, 'Type')
            % Read a link
            data = hdf5_loadLink(info);
            
        else
            error('hdf5_load:parse_argument', ...
                'The 2-nd argument must be a name of a group, dataset, or link');
        end
    end
    
elseif nargin == 3
    % Read attributes of a group or a data set, or slices of a data set
    location = varargin{1};
    try
        info = h5info(filename, location);
        if strcmp(info.Name, '/') % a special case of the root group
            info.Name = '';
        end
        
    catch ME
        if strcmp(ME.identifier, 'MATLAB:imagesci:h5info:fileOpenErr')
            ME = MException('hdf5_load:h5info', ...
                strjoin({'File', filename, 'does not exist or is not HDF5 format'}));
            
        elseif strcmp(ME.identifier, 'MATLAB:imagesci:h5info:libraryError')
            ME = MException('hdf5_load:h5info', ...
                strjoin({'H5Location', location, 'was not found in', filename, 'file'}));
        end
        
        throwAsCaller(ME);
    end
    
    if iscell(varargin{2})
        % Read slices of a data set
        slices = varargin{2};
        
        % Check if the specified location is a data set
        if ~isfield(info, 'Dataspace')
            error('hdf5_load:invalid_location', ...
                'Slicing ranges are not applicable, the location is not a data set');
        end
        
        data_size = info.Dataspace.Size;
        if length(slices) > length(data_size)
            error('hdf5_load:invalid_slicing', ...
                'A number of slicing ranges is larger than a dimention of a data set')
        end
        
        % Parse ranges
        startIndex = ones(1, length(data_size));
        nElements = Inf(1, length(data_size));
        for i = 1:length(slices)
            [startIndex(i), nElements(i)] = parse_range(slices{i}, data_size(i));
        end
        
        % Read data
        data = h5read(filename, location, startIndex, nElements);
        
    elseif any(strcmp(varargin{2}, {'-a', '-ca', '-c'}))
        % Read attributes and/or convert to matlab structures
        if any(strcmp(varargin{2}, {'-a', '-ca'}))
            load_attr = true;
        end
        if any(strcmp(varargin{2}, {'-ca', '-c'}))
            convert2matlab = true;
        end
        
        if isfield(info, 'Groups')
            % Read a group with its internal hierarchy
            [data, links] = hdf5_loadGroup(filename, info, convert2matlab, load_attr);
            data = assign_links(data, info, links, varargin{2});
            
        elseif isfield(info, 'Datatype')
            % Read a data set
            type = info.Datatype.Class;
            dset_val = hdf5_loadDataset(filename, location, type);
            if load_attr || convert2matlab
                [dset_attr, ml_class_dset] = hdf5_loadAttributes(info, convert2matlab, load_attr);
            else
                ml_class_dset = [];
            end
            
            if ~isempty(ml_class_dset)
                switch ml_class_dset
                    case 'complex'
                        dset_val = dset_val.r + 1i*dset_val.i;
                        
                    case 'cell'
                        if ~iscell(dset_val)
                            dset_val = {dset_val};
                        end
                        
                    case 'char_array'
                        dset_val = char(dset_val);
                        
                        
                    otherwise
                        conv2ml = str2func(ml_class_dset);
                        dset_val = conv2ml(dset_val);
                end
            end
            
            if load_attr
                data.Attributes = dset_attr;
                data.Value = dset_val;
                
            else
                data = dset_val;
            end
            
        elseif isfield(info, 'Type')
            % Read a link
            data = hdf5_loadLink(info, convert2matlab, load_attr);
            
        end
        
    else
        error('hdf5_load:parse_argument', ...
            'Incorrect 3-rd argument');
    end
end

function [data, links] = hdf5_loadGroup(filename, info, varargin)
import io.HDF.*

if nargin > 2
    convert2matlab = varargin{1};
    load_attr = varargin{2};
else
    convert2matlab = false;
    load_attr = false;
end

data = [];

% Collect links
links = info.Links.'; % transform to a row for easier indexing
if ~isempty(links)
    for link_ind = 1:length(links)
        links(link_ind).Name = [info.Name, '/', links(link_ind).Name];
    end
end

% Load the datasets
for dataset_ind = 1:length(info.Datasets)
    dset_info = info.Datasets(dataset_ind);
    dset_name = dset_info.Name;
    location = [info.Name, '/', dset_name];
    type = dset_info.Datatype.Class;
    
    dset_val = hdf5_loadDataset(filename, location, type);
    
    % Load attributes of a dataset
    if load_attr || convert2matlab
        [dset_attr, ml_class_dset] = hdf5_loadAttributes(dset_info, convert2matlab, load_attr);
    else
        ml_class_dset = [];
    end
    
    if ~isempty(ml_class_dset)
        switch ml_class_dset
            case 'complex'
                dset_val = dset_val.r + 1i*dset_val.i;
                
            case 'cell'
                if ~iscell(dset_val)
                    dset_val = {dset_val};
                end
                
            case 'char_array'
                dset_val = char(dset_val);
                
            otherwise
                conv2ml = str2func(ml_class_dset);
                dset_val = conv2ml(dset_val);
        end
    end
    
    if load_attr
        data.(dset_name).Attributes = dset_attr;
        data.(dset_name).Value = dset_val;
        
    else
        data.(dset_name) = dset_val;
    end
end

% Load attributes of a group
if load_attr || convert2matlab
    [group_attr, ml_class_group] = hdf5_loadAttributes(info, convert2matlab, load_attr);
    if load_attr
        data.Attributes = group_attr;
    end
else
    ml_class_group = [];
end

% Load the internal groups recursively
for group_ind = 1:length(info.Groups)
    [group_data, child_links] = hdf5_loadGroup(filename, info.Groups(group_ind), convert2matlab, load_attr);
    
    [~, group_name] = fileparts(info.Groups(group_ind).Name);
    data.(group_name) = group_data;
    
    % Aggregate links
    links = [links, child_links]; %#ok<AGROW> There shouldn't be too many links present
end

if ~isempty(ml_class_group)
    % convert the groups
    data_temp = data;
    data = [];
    if isfield(data_temp, 'Attributes')
        data.Attributes = data_temp.Attributes;
        data_temp = rmfield(data_temp, 'Attributes');
        fn = fieldnames(data_temp);
        for group_ind  = 1:length(fn)
            switch ml_class_group
                case 'cell'
                    data.Value{group_ind} = data_temp.(fn{group_ind});
                    
                case 'structure array'
                    data.Value(group_ind) = data_temp.(fn{group_ind});
                    
                otherwise
                    keyboard
            end
        end
    else
        fn = fieldnames(data_temp);
        for group_ind  = 1:length(fn)
            switch ml_class_group
                case 'cell'
                    data{group_ind} = data_temp.(fn{group_ind});
                    
                case 'structure array'
                    data(group_ind) = data_temp.(fn{group_ind});
                    
                otherwise
                    keyboard
            end
        end
    end
    
end

function [data, ml_class] = hdf5_loadAttributes(info, convert2matlab, load_attr)
data = [];
ml_class = [];
if isfield(info, 'Attributes') % info structure may not contain Attributes field
    attr_info = info.Attributes;
    for attr_ind = 1:length(attr_info)
        attr = attr_info(attr_ind);
        attr_name = attr.Name;
        if ~isvarname(attr_name)
            if ~any(strcmpi({attr_info.Name}, ['MATLAB' attr_name])) && ~strcmpi(attr_name, '_class')
                warning('Invalid attribute name! Added "MATLAB" prefix to %s.', attr_name)
                attr_name = ['MATLAB' attr_name];
            else
                error('Invalid attribute name.')
            end
        end
        if convert2matlab && strcmpi(attr_name, 'MATLAB_class')
            ml_class = attr.Value{1};
            
        elseif load_attr
            if iscell(attr.Value)
                data.(attr_name) = attr.Value{1};
            else
                data.(attr_name) = attr.Value;
            end
        end
    end
end

function data = hdf5_loadDataset(filename, location, type)
if strcmp(type, 'H5T_ENUM')
    % Workaround for a bug in h5postprocessenums (part of h5read) function
    data = read_enum(filename, location);
    
else
    data = h5read(filename, location);
    if iscell(data) && numel(data) == 1 && ischar(data{1})
        data = data{1}; % utility string unwrapping from a single cell
    end
end

function data = hdf5_loadLink(link, varargin)

if nargin > 2
    convert2matlab = varargin{1};
    load_attr = varargin{2};
else
    convert2matlab = false;
    load_attr = false;
end

switch link.Type
    case {'hard link', 'soft link'}
        filename = link.Filename;
        location = link.Value{1};
        
    case 'external link'
        filename = absolute_path(link.Value{1}, link.Filename);
        location = link.Value{2};
        
    otherwise
        error('hdf5_load:hdf5_loadLink', ...
            strjoin({'Unknown link type at', link.Name}));
end

link_info = h5info(filename, location);
if strcmp(link_info.Name, '/') % a special case of the root group
    link_info.Name = '';
end

if isfield(link_info, 'Groups')
    [data, links] = hdf5_loadGroup(filename, link_info, convert2matlab, load_attr);
    data = assign_links(data, link_info, links);
    
elseif isfield(link_info, 'Datatype')
    type = link_info.Datatype.Class;
    data = hdf5_loadDataset(filename, location, type);
    
elseif isfield(link_info, 'Type')
    data = hdf5_loadLink(link_info, convert2matlab, load_attr);
    
else
    error('hdf5_load:hdf5_loadLink', ...
        strjoin({'A link at', link_info.Name, 'must be a name of a group, dataset, or link'}));
end

function data = assign_links(data, info, links, varargin)
import io.HDF.*

if ~isempty(varargin)
    flag = varargin{1};
else
    flag = [];
end
if ~isempty(links)
    cut_start = length(info.Name) + 1;
    
    while true
        resolved_links = false(size(links));
        
        for ind = 1:length(links)
            link = links(ind);
            place = strrep(link.Name(cut_start:end), '/', '.');
            target = [];
            target_struc = [];
            if ~isempty(flag) && contains(flag, 'a') && contains(flag, 'c')
                target_struc = ['.Value'];
            end
            switch link.Type
                case {'hard link', 'soft link'}
                    try
                        parent = strsplit(link.Value{1}, '/');
                        parent = strjoin(parent(1:end-1), '/');
                        parent_info = h5info(info.Filename, parent);
                        if isfield(parent_info, 'Attributes') && ~isempty(parent_info.Attributes)
                            for ii=1:numel(parent_info.Attributes)
                                if strcmp(parent_info.Attributes(ii).Name, 'MATLAB_class') && ~isempty(flag) && contains(flag, 'c')
                                    % get pointer index
                                    pnt_indx = strsplit(link.Value{1}, '_');
                                    pnt_indx = str2double(pnt_indx(end));
                                    target_add = [];
                                    switch parent_info.Attributes(ii).Value{1}
                                        case 'cell'
                                            target_add = sprintf('{%d}', pnt_indx+1);
                                            
                                        case 'structure array'
                                            target_add = sprintf('(%d)', pnt_indx+1);
                                            
                                        otherwise
                                            keyboard
                                    end
                                    target = [strrep(parent, '/', '.') target_struc target_add];
                                    break
                                end
                            end
                        end
                        if isempty(target)
                            target = [strrep(link.Value{1}(cut_start:end), '/', '.') target_struc];
                        end
                        
                        evalc(['data', place, ' = data', target]);
                        
                    catch
                        continue % postpone this link resolution
                    end
                    
                case 'external link'
                    ext_link = absolute_path(link.Value{1}, info.Filename);
                    
                    % make sure to reference the same variable in evalc!
                    if ~isempty(flag)
                        target_data = hdf5_load(ext_link, link.Value{2}, flag); %#ok<NASGU>
                    else
                        target_data = hdf5_load(ext_link, link.Value{2});
                    end
                    evalc(['data', place, ' = target_data']);
                    
                otherwise
                    error('hdf5_load:assign_links', ...
                        strjoin({'Unknown link type at', place}));
            end
            
            resolved_links(ind) = true;
        end
        
        if all(resolved_links)
            % all links have been assigned
            return
        end
        
        if ~any(resolved_links)
            % none of the links has been assigned in this iteration
            error('hdf5_load:assign_links', ...
                strjoin({'Cannot assign link(s) at', ''}));
        end
        
        links = links(~resolved_links);
    end
end

function filepath = absolute_path(filepath, current_filepath)
if ~startsWith(filepath, '/')
    path = fileparts(current_filepath);
    filepath = fullfile(path, filepath);
end

function [startVal, nVals] = parse_range(valRange, maxVal)
if isempty(valRange) % empty
    startVal = 1;
    nVals = Inf;
    
elseif isscalar(valRange) % single value
    if valRange <= -1
        valRange = maxVal + valRange + 1;
    end
    startVal = valRange;
    nVals = 1;
    
elseif isvector(valRange) && numel(valRange) == 2 % vector with two values
    if valRange(1) <= -1
        if isinf(valRange(1))
            valRange(1) = 1; % = -Inf
        else
            valRange(1) = maxVal + valRange(1) + 1;
        end
    end
    startVal = valRange(1);
    
    if valRange(2) <= -1
        if isinf(valRange(2))
            valRange(2) = 1; % = -Inf
        else
            valRange(2) = maxVal + valRange(2) + 1;
        end
    end
    nVals = valRange(2) - startVal + 1;
    
else
    error('hdf5_load:parse_range', ...
        'A range should be specified with <= 2 parameters');
end

if startVal < 1 || startVal > maxVal || nVals < 1 || (nVals > maxVal && ~isinf(nVals))
    error('hdf5_load:parse_range', ...
        'The resulting range is out of data borders');
end

function data = read_enum(filename, location)
file_id = H5F.open(filename);
dset_id = H5D.open(file_id, location);
type_id = H5D.get_type(dset_id);

data = H5D.read(dset_id); % numerical member of enumeration
data = H5T.enum_nameof(type_id, data); % associated symbol name

H5T.close(type_id);
H5D.close(dset_id);
H5F.close(file_id);

