% FIND_PROJECTION_FILES_NAMES_APS find names of ptychography projections given
% the tomo param structure, if more files are found, return path to the
% newest one. Written by YJ for ML reconstructions
%
% filename = find_projection_files_names(par, scan_num)
% 
% Inputs
%   **par                tomo par structure 
%   **scan_num           scan number 
% *returns*
%   ++filename           filename of the found scan 

function filename = find_projection_files_names_aps(par, scan_num)
    %bug?
    %%%data_in_subfolders = exist(fullfile(par.analysis_path, 'S00000-00999'), 'dir');
    data_in_subfolders = exist(fullfile(par.analysis_path, 'S02000-02999'), 'dir');
    
    if data_in_subfolders
        % new x12sa path format
        % bug?
        %%%path = fullfile(par.analysis_path, utils.compile_x12sa_dirname(scan_num) , [par.fileprefix '*' par.filesuffix '*.' par.file_extension]); 
    	%modified by YJ
        path = fullfile(par.analysis_path, utils.compile_x12sa_dirname(scan_num) , [sprintf('S%05i',scan_num), par.fileprefix '*' par.filesuffix '*.' par.file_extension]); 
        
    else
        % compatibility option for original analysis paths 
        
        %path = fullfile(par.analysis_path,sprintf(par.scan_string_format,scan_num),[par.fileprefix '*' par.filesuffix '*.' par.file_extension]); 
        
       
        
    end
    
    while contains(path, '**')
        path = replace(path, '**', '*'); % prevent failure when path string contains multiple asterix
    end
    path = utils.abspath(path); 
    path = dir(path); 
    if isempty(path)
        filename = []; 
        return
    end
    % take the last file fitting the constraints 
    [~,ind]=sort([path.datenum]);
    filename = fullfile(path(ind(1)).folder, path(ind(1)).name);
end