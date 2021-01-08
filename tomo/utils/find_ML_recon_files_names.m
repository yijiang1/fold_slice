% FIND_ML_RECON_FILES_NAMES find names of ML ptychographic reconstructions given
% the param structure defined by YJ for the APS datasets.
% user can search for differernt # of probes and OPR modes
%
% filename = find_projection_files_names(par, scan_num)
% 
% Inputs
%   **par                tomo par structure 
%   **scan_num           scan number 
% *returns*
%   ++filename           filename of the found scan 


function [filename,method,roi,scanNo] = find_ML_recon_files_names(par, scan_num)
    N_roi = linspace(1,length(par.MLrecon.roi),length(par.MLrecon.roi));
    [N_roi_temp, Nprobe_s_temp,var_probe_modes_s_temp] = ndgrid(N_roi,par.MLrecon.Nprobes, par.MLrecon.var_probe_modes);
    N_roi_temp = reshape(N_roi_temp,[],1);
    Nprobe_s_temp = reshape(Nprobe_s_temp,[],1);
    var_probe_modes_s_temp = reshape(var_probe_modes_s_temp,[],1);
    filename = [];
    scanNo = num2str(scan_num);
    for i=1:length(Nprobe_s_temp)
        roi = par.MLrecon.roi{N_roi_temp(i)};
        Nprobe = Nprobe_s_temp(i);
        vp_modes = var_probe_modes_s_temp(i);
        %method = sprintf(par.MLrecon.method,Nprobe,vp_modes,Nprobe,vp_modes);
        method = '';
        for j=1:length(par.MLrecon.method)
            method = strcat(method,sprintf(par.MLrecon.method{j},Nprobe,vp_modes));
        end
        
        filename_temp = strcat(par.MLrecon.path, sprintf(par.scan_string_format, scan_num),roi, method,'/Niter',num2str(par.MLrecon.Niter),'.mat');
        %disp('filename_temp:')
        %disp(filename_temp)
        if ~isempty(dir(filename_temp))
            d = dir(filename_temp);
            %filename = filename_temp;
            filename = strcat(d(1).folder,'/',d(1).name);
            %disp('filename:')
            %disp(filename)
            %disp(strcat(d.folder,'/',d.name))

            %disp(strcat(sprintf(par.scan_string_format, scan_num),'--',method,'Niter',num2str(par.MLrecon.Niter),'.mat'))
            break
        end
    end
    
    %{
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
        %path = fullfile(par.analysis_path,sprintf('S%05i',scan_num),[par.fileprefix '*' par.filesuffix '*.' par.file_extension]); 
        % path = strcat(par.analysis_path,sprintf(par.scan_string_format,scan_num),'.', par.file_extension); 
        %generate a list of possibe recon parameters
        %current support: number of probe modes; number of OPR modes
        [Nprobe_s_temp,var_probe_modes_s_temp] = meshgrid(par.MLrecon.Nprobes, par.MLrecon.var_probe_modes);
        Nprobe_s_temp = reshape(Nprobe_s_temp,[],1);
        var_probe_modes_s_temp = reshape(var_probe_modes_s_temp,[],1);
        for i=1:length(Nprobe_s_temp)
            Nprobe = Nprobe_s_temp(i);
            vp_modes = var_probe_modes_s_temp(i);
            method = sprintf(par.MLrecon.method,Nprobe,vp_modes,Nprobe,vp_modes);
            filename = strcat(par.MLrecon.path, sprintf(par.scan_string_format, scan_num),par.MLrecon.roi, method,'Niter',num2str(par.MLrecon.Niter),'.mat');
            if ~isempty(dir(filename))
                disp(strcat(sprintf(par.scan_string_format, scan_num),'--',method,'Niter',num2str(par.MLrecon.Niter),'.mat'))
                break
            end
        end
    end
    %}
    %{
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
    %}
end