clear variables
addpath(strcat(pwd,'/utils/'))
addpath(core.find_base_package)

%%
%scanNo_s = [276]; 

gpu_id = 2; %gpu0
%gpu_id = 3; %gpu1
%gpu_id = 5; %gpu3

%gpu_id = 1; %gpu4
%gpu_id = 6; %gpu5
%gpu_id = 7; %gpu6
%gpu_id = 8; %gpu7
%gpu_id = 9; %gpu8
%gpu_id = 10; %gpu9

scanNo_s = [94:252];
executed_template = 'ptycho_bnp_19c3_Lin_single';
base_path = '/mnt/micdata1/bnp/2019-3/Lin/results/Lin_ML2/';
dataName = 'data_roi0_Ndp64_dp.hdf5'; %skip recon if data doesn't exist
N_attampt = 1; %allow multiple attempts
log_title = 'BNP 19C3 Lin';
log_notes = ''; %write extra info

%%
%{
try
    run(executed_template)
catch
    error(['Fail to run recon template: ',executed_template '.m'])
end
%}
reconLog = strcat(base_path,'ML_recon_log/');
if ~exist(reconLog,'dir'); mkdir(reconLog); end

while true
    all_recon_finished = true;
    for i=1:length(scanNo_s)
        scanNo = scanNo_s(i);

        %check if data is processed already
        reconLogName = strcat(reconLog,'scan',num2str(scanNo_s(i)),'*.txt');
        if ~isempty(dir(reconLogName)); continue; end
        %check if data exist
        dataPath = strcat(base_path,'scan',num2str(scanNo_s(i)),'/',dataName);
        %disp(dataPath)
        if ~exist(dataPath,'file'); continue; end
        %disp(scanNo_s(i))
        
        %% write new recon log
        ct = clock;
        reconLogName = strcat(reconLog,'scan',num2str(scanNo_s(i)),'.txt');
        fileID = fopen(reconLogName,'a');
        fprintf(fileID,[log_title,'\n']);
        if ~isempty(log_notes); fprintf(fileID,[log_notes,'\n']); end
        fprintf(fileID,['User: ',io.get_user_name(),'\n']);
        fprintf(fileID,['Server: ',io.get_host_name(),'\n']);
        fprintf(fileID,'GPU_id=%d\n',gpu_id);
        fprintf(fileID,'Start time: %d-%d-%d %.2d:%.2d:%.2d\n',ct(2),ct(3),ct(1),ct(4),ct(5),round(ct(6)));
        fclose(fileID);
        all_recon_finished = false;
        try
            reset(gpuDevice( gpu_id ))
        catch
        end
        %begin recon
        for ia = 1:N_attampt
            try
                run(executed_template) %initialize p and eng
                % Run the reconstruction
                tic
                out = core.ptycho_recons(p);
                toc
                if out.recon_success
                    recon_success = true;
                    break
                end
            catch
                recon_success = false;
            end
        end
        %{
        if recon_success
            reconLogName_done = strcat(reconLog,'scan',num2str(scanNo_s(i)),'_done.txt');
            movefile(reconLogName_recon, reconLogName_done)
        else
            reconLogName_done = strcat(reconLog,'scan',num2str(scanNo_s(i)),'_fail.txt');
            movefile(reconLogName_recon, reconLogName_done)
        end
        %}
        ct = clock;
        fileID = fopen(reconLogName,'a');
        fprintf(fileID,'End time: %d-%d-%d %.2d:%.2d:%.2d\n',ct(2),ct(3),ct(1),ct(4),ct(5),round(ct(6)));
        if recon_success
             fprintf(fileID,['Recon status: ','success\n']);
        else
             fprintf(fileID,['Recon status: ','fail\n']);
        end
        fclose(fileID);
    end
    %pause(60)
    if all_recon_finished
        break
    end
end

%% send email notification
try
    p.io.script_name = mfilename;
    io.email_notification(p)
catch
end
%clear variables
