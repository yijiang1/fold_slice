% FIND_PROJECTION_FILES_NAMES find names of ptychography projections given
% the tomo param structure, if more files are found, return path to the
% newest one 
%
% filename = find_projection_files_names(par, scan_num)
% 
% Inputs
%   **par                tomo par structure 
%   **scan_num           scan number 
% *returns*
%   ++filename           filename of the found scan 


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


function filename = find_projection_files_names(par, scan_num)
    %bug?
    data_in_subfolders = exist(fullfile(par.analysis_path, 'S00000-00999'), 'dir');
    %data_in_subfolders = exist(fullfile(par.analysis_path, 'S02000-02999'), 'dir');
    %disp(data_in_subfolders)
    if data_in_subfolders
        % new x12sa path format
        % bug?
        path = fullfile(par.analysis_path, utils.compile_x12sa_dirname(scan_num) , [par.fileprefix '*' par.filesuffix '*.' par.file_extension]); 
    	
        %modified by YJ
        %path = fullfile(par.analysis_path, utils.compile_x12sa_dirname(scan_num) , [sprintf('S%05i',scan_num), par.fileprefix '*' par.filesuffix '*.' par.file_extension]); 
    else
        % compatibility option for original analysis paths 
        path = fullfile(par.analysis_path,sprintf('S%05i',scan_num),[par.fileprefix '*' par.filesuffix '*.' par.file_extension]); 
    end
    %disp(path)

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