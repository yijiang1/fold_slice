%LOAD_PTYCHO_RECONS Load data from cxs/h5 or mat file and return it as 
% structure, single dataset or directly into the workspace.
% An additional argument can be passed to select subsections of the data.
% Loading single datasets is only supported for at least 2 output
% arguments.
%
%   file...     path to cxs/h5 or mat file
%
%   *optional*
%   section...  'full', 'probe', 'object', 'recon' or 'p' to select 
%               subsections of the data; default: 'full'
%
%   EXAMPLES:
%     %% recommended usage %%
%       % load into a structure
%       S = load_ptycho_recons('./recon.h5');
%
%       % load a subset
%       S = load_ptycho_recons('./recon.h5', 'probe');
%
%       % load into single datasets
%       [object, probe, p] = load_ptycho_recons('./recon.h5');
%
%     %% not recommended, only works in 'base' workspace %%
%       % load directly into workspace
%       load_ptycho_recons('./recon.h5');
%
%
%   full = object, probe (current scan) and p
%   recon = object and probe (current scan)
%   probe = probe (current scan)
%   object = object (current scan)
%

function varargout  = load_ptycho_recons( filename_with_path, varargin )

import io.HDF.hdf5_load

varargout = {};

if ~ischar(filename_with_path)
   error('First argument has to be string') 
end

filename_with_path = utils.abspath(filename_with_path); 

if ~exist(filename_with_path, 'file')
    error('Could not find reconstruction file %s', filename_with_path)
end

if nargin > 1
    switch varargin{1}
        case {'pr'; 'probe'; 'probes'}
            section = 'probe';
        case {'ob'; 'obj'; 'objects'}
            section = 'object';
        otherwise
            section = varargin{1};
    end
else
    section = 'full';
end

if ~nargout
    output = 0;
elseif nargout >=2
    output = 2;
else
    output = 1;
end

    function assign_struct(val, val_name)
        switch output
            case 1
                varargout{1}.(val_name) = val;
            case 2
                varargout{end+1} = val;
            otherwise 
                assignin('base', val_name, val);
        end
    end

    function assign_val(struc)
        switch output
            case 1
                varargout{1} = struc;
                
            case 2
                if isfield(struc, 'object')
                    varargout{end+1} = struc.object;
                end
                if isfield(struc, 'probe')
                    varargout{end+1} = struc.probe;
                end
                if isfield(struc, 'p')
                    varargout{end+1} = struc.p;
                end
                    
            otherwise
                fn = fieldnames(struc);
                for ii=1:length(fn)
                    assignin('base', fn{ii}, struc.(fn{ii}))
                end
        end
    end


% check if it is a .mat file or a .cxs file
[~, ~, ext] = fileparts(filename_with_path);
switch ext
    case '.mat'
        switch section
            case 'recon'
                S = load(filename_with_path, 'object', 'probe');
                assign_val(S);

            case 'full'
                S = load(filename_with_path);
                assign_val(S);

            case 'object'
                S = load(filename_with_path, 'object');
                assign_val(S);

            case 'probe'
                S = load(filename_with_path, 'probe');
                assign_val(S);

            case 'p'
                S = load(filename_with_path, 'p');
                assign_val(S);

            otherwise
                error('Unknown data section %s', section);
        end
        
    case {'.cxs','.h5'}
        if io.HDF.hdf5_dset_exists(filename_with_path, 'object', '/reconstruction', true)
            h5_path = '/reconstruction';
        else
            h5_path = '';
        end
        
        % reconstruction
        switch section
            case 'recon'
                % load object
                h = hdf5_load(filename_with_path, [h5_path '/object']);
                assign_struct(load_data_cell(h), 'object');
                % load probe
                h = hdf5_load(filename_with_path, [h5_path '/probes']);
                assign_struct(load_data_cell(h), 'probe');
            case 'full'
                % load object
                h = hdf5_load(filename_with_path, [h5_path '/object']);
                assign_struct(load_data_cell(h), 'object');
                % load probe
                h = hdf5_load(filename_with_path, [h5_path '/probes']);
                assign_struct(load_data_cell(h), 'probe');
                % load p
                p = convert2p(hdf5_load(filename_with_path, '/reconstruction/p', '-c'));
                if io.HDF.hdf5_dset_exists(filename_with_path, 'meta_all', '/measurement', true)
                    p.meta = hdf5_load(filename_with_path, '/measurement/meta_all', '-c');
                elseif io.HDF.hdf5_dset_exists(filename_with_path, 'spec_all', '/measurement', true)
                    p.meta = hdf5_load(filename_with_path, '/measurement/spec_all', '-c');
                end
                assign_struct(p, 'p');                
            case 'object'
                % load object
                h = hdf5_load(filename_with_path, [h5_path '/object']);
                assign_struct(load_data_cell(h), 'object');
            case 'probe'
                % load probe
                h = hdf5_load(filename_with_path, [h5_path '/probes']);
                assign_struct(load_data_cell(h), 'probe');
            case 'p'
                % load p
                p = convert2p(hdf5_load(filename_with_path, '/reconstruction/p', '-c'));
                if io.HDF.hdf5_dset_exists(filename_with_path, 'meta_all', '/measurement', true)
                    p.meta = hdf5_load(filename_with_path, '/measurement/meta_all', '-c');
                elseif io.HDF.hdf5_dset_exists(filename_with_path, 'spec_all', '/measurement', true)
                    p.meta = hdf5_load(filename_with_path, '/measurement/spec_all', '-c');
                end
                assign_struct(p, 'p');
            otherwise
                error('Unknown data section %s', section);
        end
                
        
        
    otherwise
        error('Unknown ptycho datatype %s.', ext)
end



end

function tmp = load_data_cell(h)

    fn = fieldnames(h);
    num_end = str2double(subsref(strsplit(fn{1}, '_'), struct('type', '{}', 'subs',{{length(strsplit(fn{1},'_'))}})));
    if length(fn)==2 && (strcmpi(fn{1}, 'i') || strcmpi(fn{1}, 'r'))
        tmp = permute(h.r + 1i*h.i, [2,1,3,4]);
    elseif isnumeric(num_end) && ~isnan(num_end)
        for ii=1:length(fn)
            if isstruct(h.(fn{ii}))
                tmp{ii} = load_data_cell(h.(fn{ii}));
            else
                if isnumeric(h.(fn{ii}))
                    tmp{ii} = double(h.(fn{ii}));
                else
                    tmp{ii} = h.(fn{ii});
                end
            end
        end
%         tmp = h;
    else
        for ii=1:length(fn)
            if isstruct(h.(fn{ii}))
                tmp.(fn{ii}) = load_data_cell(h.(fn{ii}));
            else
                if isnumeric(h.(fn{ii}))
                    tmp.(fn{ii}) = double(h.(fn{ii}));
                else
                    tmp.(fn{ii}) = h.(fn{ii});
                end
            end
        end
    end

end

function tmp = convert2p(h)
    
    fn = fieldnames(h);
    for ii=1:length(fn)
        if isstruct(h.(fn{ii}))
            h.(fn{ii}) = load_data_cell(h.(fn{ii}));
        elseif isnumeric(h.(fn{ii}))
            h.(fn{ii}) = double(h.(fn{ii}));                
        else
            continue;
        end
    end
    tmp = h;
    
    % object
    for ii=1:length(h.objects)
        tmp.object{ii} = permute(h.objects{ii}, [2 1 3 4]);
    end
    tmp = rmfield(tmp, 'objects');
    
    % probes
    pr = tmp.probes;
    tmp.probes = [];
    for ii=1:length(pr)
        tmp.probes(:,:,ii,:) = permute(pr{ii}, [2 1 3 4]);
    end
    
    % positions
    tmp.positions = transpose(tmp.positions);
    tmp.positions_real = transpose(tmp.positions_real);
    tmp.positions_orig = transpose(tmp.positions_orig);
    
    % ctr
    tmp.ctr = transpose(tmp.ctr);
            
end

