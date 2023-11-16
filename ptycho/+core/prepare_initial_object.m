%PREPARE_INITIAL_OBJECT 
% prepare an initial guess for the object reconstruction
%
% ** p      p structure
%
% returns:
% ++ p      p structure
%
% see also: core.prepare_initial_guess
%

function [ p ] = prepare_initial_object( p )
import utils.verbose 
import utils.crop_pad 
import utils.interpolateFT 

if p.model_object
    switch p.model.object_type
        case 'rand'
            verbose(2,'Using random object as initial guess.')
        case 'amplitude'
            % create dummy object which will be overwritten once the
            % prepared data is available
            assert(p.fourier_ptycho, 'An initial guess based on the prepared data is only suited for Fourier ptychography.')
    end
    for obnum = 1:p.numobjs
        p.object{obnum} = (1+1i*1e-6*rand([p.object_size(obnum,:) p.object_modes])).*ones([p.object_size(obnum,:) p.object_modes]);
    end
else
    if isfield(p, 'initial_iterate_object')
        warning('Loading initial object guess from file given by p.initial_iterate_object_file.')
    end
    verbose(2,'Using loaded object as initial guess.')

    if numel(p.initial_iterate_object_file) ~= p.numobjs
        verbose(2,'Number of initial iterate files and number of objects does not match')
        for ii=numel(p.initial_iterate_object_file):p.numobjs
            p.initial_iterate_object_file{ii} = p.initial_iterate_object_file{end};
        end
    end
    
    % make a bit smarter the use of initial_iterate_object_file and allow
    % some automatic patten filling + file search 
    for obnum = unique(p.share_object_ID)
        % if string allows it, fill in the scan numbers 
        p.initial_iterate_object_file{obnum} = sprintf(p.initial_iterate_object_file{obnum}, p.scan_number(obnum)); 
        if contains(p.initial_iterate_object_file{obnum}, '*') % if string contains wild character *, try to find the file 
           fpath = dir(p.initial_iterate_object_file{obnum}) ; 
           if isempty(fpath)
               warning('No file corresponding to pattern %s was found, using random initial guess', p.initial_iterate_object_file{obnum})
               p.object{obnum} = (1+1i*1e-6*rand([p.object_size(obnum,:) p.object_modes])).*ones([p.object_size(obnum,:) p.object_modes]);
               p.initial_iterate_object_file{obnum} = []; 
               continue
           elseif length(fpath) > 1
               warning('Too many files corresponding to pattern %s were found, using the last', p.initial_iterate_object_file{obnum})
               fpath = fpath(end); 
           end
           p.initial_iterate_object_file{obnum} = [fpath.folder,'/',fpath.name]; 
        end
    end
    
    % load data from disk
    for ii = unique(p.share_object_ID) % avoid loading datasets twice
        if isempty(p.initial_iterate_object_file{ii})
            continue
        end
        if ~exist(p.initial_iterate_object_file{ii}, 'file')
            error(['Did not find initial iterate: ' p.initial_iterate_object_file{ii}])
        end
       
        verbose(2,'Loading object %d from: %s',ii,p.initial_iterate_object_file{ii})
        S = io.load_ptycho_recons(p.initial_iterate_object_file{ii});
        object = double(S.object);
        % reinterpolate to the right pixel size 
        if isfield(S, 'p') && any(S.p.dx_spec ~= p.dx_spec)
            verbose(2, 'Warning: Reinterpolate loaded object to new pixels size')
            object = interpolateFT(object, ceil(size(object(:,:,1)).*S.p.dx_spec./p.dx_spec));
        end

        %%% check the object size
        if ~isequal(size(squeeze(object(:,:,1))), squeeze(p.object_size(ii,:)))
            % if the loaded dataset does not have the expected object size,
            % crop/pad it to p.object_size
            verbose(2, 'Warning: Object taken from file %s does not have the expected size of %d x %d.', ...
            p.initial_iterate_object_file{ii}, p.object_size(ii,1), ...
            p.object_size(ii,2))
            p.object{ii} = crop_pad(object, p.object_size(ii,:));
        else
            % the the object sizes are the same, just copy everything
            % to p.object
            p.object{ii} = object;
        end
        
        % now let's check the object modes
        mode_diff = p.object_modes-size(object,3);
        if mode_diff > 0
            % add (random) object modes
            p.object{ii}(:,:,size(object,3)+1:p.object_modes,:) = (1+1i*1e-6*rand([p.object_size(ii,:) mode_diff])).*ones([p.object_size(ii,:) mode_diff]);
        elseif mode_diff < 0
            % modified by YJ. keep all layers for multi-layer object
            if isfield(p,'multiple_layers_obj') && p.multiple_layers_obj
                if isfield(p,'sum_obj_layers') && p.sum_obj_layers
                    verbose(2, 'Sum all layers in the initial object.')
                    object_temp = prod(p.object{ii},3);
                else
                    %add an extra axis that is needed by GPU_MS
                    object_temp(:,:,1,:) = p.object{ii}; 
                end
                p.object{ii} = object_temp;
            else
                % remove object modes
                p.object{ii}(:,:,p.object_modes+1:size(object,3),:) = [];
            end
        end
    end
end

end
