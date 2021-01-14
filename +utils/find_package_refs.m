%FIND_PACKAGE_REFS update references in <path> and its subfolders to
%package structure in <base_path>. 
%   
%   base_path...        repository with new package structure
%   path...             repository which needs to be updated
%   
%   *optional*          given as name/value pair
%   extension...        file extension; default '.m'
%   recursive...        recursive behavior; default false
%   show_files...       show file names, otherwise progressbar; default false
%   filename...         change output file name and path; default
%                       ./references.txt
%   
%   Example: 
%       find_package_refs('./cSAXS_matlab_base', './cSAXS_matlab_ptycho')
%       find_package_refs('./cSAXS_matlab_base', './cSAXS_matlab_ptycho', 'recursive', false);
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


function find_package_refs(base_path, path, varargin )   

% check path
if ~exist(path)
    error('Could not find %s', path)
end
if ~exist(base_path)
    error('Could not find %s', base_path)
end

% set default values
extension = '.m';
recursive = false;
show_files = true;
filename_with_path = '../references.txt';

% parse the variable input arguments vararg = cell(0,0);
if ~isempty(varargin)
    for ind = 1:2:length(varargin)
        name = varargin{ind};
        value = varargin{ind+1};
        switch lower(name)
            case 'extension'
                extension = value;
            case 'recursive'
                recursive = value;
            case 'show_files'
                show_files = value;
            case 'filename'
                filename_with_path = value;
        end
    end
end

% avoid overwriting files
if exist(filename_with_path,'file')
    disp(['File ' filename_with_path ' exists,' ])
    userans = input(['Do you want to overwrite (y/N)? '],'s');
    if strcmpi(userans,'y')
        disp(['Saving to  ' filename_with_path]);
    else
        disp(['Did not save ' filename_with_path])
        return
    end
else
    display(['Saving to  ' filename_with_path]);
end


fileID = fopen(filename_with_path,'w');

file_list = [];




% get the target file list
if recursive
    [~, target_fl] = system(['find ' path ' -name "*' extension '"']);
    target_fl = strsplit(target_fl, '\n');
else
    [target_fl_temp] = dir([path '/*' extension]);
    for ii=1:length(target_fl_temp)
        target_fl{ii} = [path '/' target_fl_temp(ii).name];
    end
    
end

% get the package file list
[~, base_fl] = system(['find ' base_path ' -name "*' extension '"']);
base_fl = strsplit(base_fl, '\n');
for ii=1:length(base_fl)
    pckg_name = {};
    substr = strsplit(base_fl{ii}, '/');
    fn = substr{end};
    if isempty(fn)
        continue
    else
        % get the updated package name
        for jj=1:length(substr)
            try
                if strcmp(substr{jj}(1), '+')
                    pckg_name{end+1} = substr{jj}(2:end);
                    
                end
            catch
                continue
            end
        end
    end
    
    pckg_name_full = strjoin(pckg_name, '.');

    if show_files
        fprintf('-- Updating references to file %s.\n', fn);
    end
    
    % call external function and update reference
    for kk=1:length(target_fl)
        if ~isempty(target_fl{kk})
            temp_fn = strsplit(target_fl{kk}, '/');
            temp_fn = temp_fn{end};
            if ~isfield(file_list, temp_fn(1:end-length(extension)))
                file_list.(temp_fn(1:end-length(extension))).pckgs = [];
                file_list.(temp_fn(1:end-length(extension))).files = [];
                file_list.(temp_fn(1:end-length(extension))).subfunctions = [];
            end
            
            [~, cnt] = system(['grep -n ' fn(1:end-length(extension)) ' ' target_fl{kk} '| wc -l']);
            count = str2double(cnt);
            
            if strcmp(fn, temp_fn) && count<=1
                fprintf('Skipping %s\n', fn)
                continue
            elseif count >=1
                [~, cnt] = system(['grep -n ''^function'' ' target_fl{kk} '| wc -l']);
                count = str2double(cnt) -1;
                file_list.(temp_fn(1:end-length(extension))).pckgs{end+1} = pckg_name_full;
                file_list.(temp_fn(1:end-length(extension))).files{end+1} = [pckg_name_full '.' fn(1:end-length(extension))];
                file_list.(temp_fn(1:end-length(extension))).subfunctions = (count>0)*count;
            end
        end
    end
    if ~show_files
        utils.progressbar(ii, length(base_fl))
    end
end

fn = fieldnames(file_list);
for ii=1:length(fn)
    if ~isempty(unique(file_list.(fn{ii}).pckgs(:)))
        fprintf(fileID, [fn{ii} '\n']);
        fprintf(fileID, [strjoin(unique(file_list.(fn{ii}).pckgs(:)), '\t') '\n']);
        fprintf(fileID, [strjoin(unique(file_list.(fn{ii}).files(:)), '\t') '\n']);
        fprintf(fileID, ['Subfunctions: ' num2str(file_list.(fn{ii}).subfunctions)]);
        fprintf(fileID, ['\n\n']);

    end
end


fclose(fileID);
    

end

