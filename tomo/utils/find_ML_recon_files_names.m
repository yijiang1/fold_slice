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


function [filename, method, roi, scanNo] = find_ML_recon_files_names(par, scan_num)
    N_roi = linspace(1,length(par.MLrecon.roi),length(par.MLrecon.roi));
    filename = [];
    scanNo = num2str(scan_num);

    if isempty(par.MLrecon.var_probe_modes) % no variable probe correction
        [N_roi_temp, Nprobe_s_temp] = ndgrid(N_roi, par.MLrecon.Nprobes);
        N_roi_temp = reshape(N_roi_temp,[],1);
        Nprobe_s_temp = reshape(Nprobe_s_temp,[],1);
        for i=1:length(Nprobe_s_temp)
            roi = par.MLrecon.roi{N_roi_temp(i)};
            Nprobe = Nprobe_s_temp(i);
            method = '';
            for j=1:length(par.MLrecon.method)
                method = strcat(method,sprintf(par.MLrecon.method{j},Nprobe));
            end
            filename_temp = strcat(par.MLrecon.path, sprintf(par.scan_string_format, scan_num),roi, method,'/Niter',num2str(par.MLrecon.Niter),'.mat');
            if ~isempty(dir(filename_temp))
                d = dir(filename_temp);
                filename = strcat(d(end).folder,'/',d(end).name);
                break
            end
        end
    else
        [N_roi_temp, Nprobe_s_temp,var_probe_modes_s_temp] = ndgrid(N_roi, par.MLrecon.Nprobes, par.MLrecon.var_probe_modes);
        N_roi_temp = reshape(N_roi_temp,[],1);
        Nprobe_s_temp = reshape(Nprobe_s_temp,[],1);
        var_probe_modes_s_temp = reshape(var_probe_modes_s_temp,[],1);
        
        for i=1:length(Nprobe_s_temp)
            roi = par.MLrecon.roi{N_roi_temp(i)};
            Nprobe = Nprobe_s_temp(i);
            vp_modes = var_probe_modes_s_temp(i);
            method = '';
            for j=1:length(par.MLrecon.method)
                method = strcat(method,sprintf(par.MLrecon.method{j},Nprobe,vp_modes));
            end
            filename_temp = strcat(par.MLrecon.path, sprintf(par.scan_string_format, scan_num),roi, method,'/Niter',num2str(par.MLrecon.Niter),'.mat');
            if ~isempty(dir(filename_temp))
                d = dir(filename_temp);
                filename = strcat(d(end).folder,'/',d(end).name);
                break
            end
        end
    end
    
    if isfield(par.MLrecon, 'print_file_path') && par.MLrecon.print_file_path
        disp(filename)
    end
end
