%PTYCHO_EXIT 
%onCleanup function for core.ptycho_recons
% 
% ** p          p structure
%
% see also: core.ptycho_recons
% modified by YJ

function ptycho_exit(p)
import utils.*

if ~isfield(p.io,'data_descriptor')
    p.io.data_descriptor = '';
end

if isfield(p.getReport, 'crashed') && p.getReport.crashed
    fprintf('\n\n\n###############################################################\n')
    fprintf('########################## PONG! ##############################\n')
    fprintf('###############################################################\n')
    disp([p.getReport.ME.getReport '\n\n\n']);

    fprintf('Reconstruction stopped! Set verbose level to >3 for debugging. \n')
    if ~p.getReport.completed
        if ~isempty(p.io.phone_number) && p.io.send_crashed_recon_SMS
            %modified by YJ
        
            subject = 'Ptycho recon crashed!';
            message = sprintf('User %s''s ptycho recon of %s scan%s crashed on %s', io.get_user_name(), ...
                p.io.data_descriptor, num2str(p.scan_number), io.get_host_name());
            if isfield(p.engines{1},'use_gpu') && p.engines{1}.use_gpu 
                message = strcat(message,'(GPU id ',num2str(p.engines{1}.gpu_id),')');
            end
            io.sendSMS(p.io.phone_number, subject, message);
        end
        if isfield(p.queue, 'file_this_recons')
            disp('test test test test')
            try
                fprintf('Moving queue file file back to %s.\n', fullfile(p.queue.path, p.queue.file_this_recons))
                
                % get log file name
                [~, ~, fext] = fileparts(p.queue.file_this_recons);
                log_dir = fullfile(p.queue.path, 'failed', 'log');
                if ~exist(log_dir, 'dir')
                    mkdir(log_dir)
                end
                log_file = fullfile(log_dir, strrep(p.queue.file_this_recons, fext, '.log'));
                
                % check if log file extists and update its content; move
                % queue file back to in_progess 
                if exist(log_file, 'file')
                    fid = fopen(log_file);
                    log_line = fgetl(fid);
                    fclose(fid);
                    log_int = strtrim(strsplit(log_line, ':'));
                    log_int = log_int{end};
                    log_int = str2double(log_int);
                    if ~isfield(p, 'queue_max_attempts')
                        p.queue.max_attempts = 5;
                        fprintf('Code crashed before parsing p.queue.max_attempts.\n')
                    end
                    if log_int >= p.queue.max_attempts
                        io.movefile_fast(fullfile(p.queue.path,'in_progress', p.queue.file_this_recons),fullfile(p.queue.path, 'failed', p.queue.file_this_recons))
                        fid = fopen(log_file, 'w');
                        fprintf(fid, [p.getReport.ME.getReport '\n\n\n']);
                        fprintf('Failed more than %u times. Moving file to ''failed''.\n', p.queue.max_attempts);
                        fclose(fid);
                        if ~isempty(p.io.phone_number) && p.io.send_failed_scans_SMS
                            io.sendSMS(p.io.phone_number, sprintf('Failed to reconstruct scan %s. I will move it to "failed".', num2str(p.scan_number)), 'sleep', p.SMS_sleep, 'logfile', fullfile(log_dir, 'sendSMS.log'));
                        end
                    else    
                        io.movefile_fast(fullfile(p.queue.path,'in_progress', p.queue.file_this_recons),fullfile(p.queue.path, p.queue.file_this_recons));
                        fid = fopen(log_file, 'w');
                        fprintf(fid, 'failed attempts: %u\n\n', log_int+1);
                        fprintf(fid, [p.getReport.ME.getReport '\n\n\n']);
                        fclose(fid);
                    end
                else
                    fid = fopen(log_file, 'w');
                    fprintf(fid, 'failed attempts: 1');
                    fclose(fid);
                    io.movefile_fast(fullfile(p.queue.path,'in_progress', p.queue.file_this_recons),fullfile(p.queue.path, p.queue.file_this_recons))
                end
                    
            catch
                fprintf('Failed to move file back to queue search path.\n')
            end
        end
        if isfield(p.queue, 'lockfile') 
            if isempty(p.queue.lockfile)
                if verbose > 2
                    p.queue.lockfile = false;
                else
                    p.queue.lockfile = true;
                end
            end
            if p.queue.lockfile
                if isempty(p.save_path{1})
                    try
                        for ii = 1:length(p.scan_number)
                            p.scan_str{ii} = sprintf(p.scan_string_format, p.scan_number(ii));        % Scan string
                        end
                        p = core.ptycho_prepare_paths(p);
                    catch
                        fprintf('Could not find lock file. \n')
                    end
                end
                for ii=1:length(p.save_path)
                    lock_filename = [p.save_path{ii} '/' p.run_name '_lock'];
                    if exist(lock_filename, 'file')
                        try
                            unix(['rm ' lock_filename]);
                            fprintf('Removing lock file %s\n',lock_filename)
                        catch
                            fprintf('Removing lock file %s failed\n',lock_filename)
                        end
                    end
                end
            end
        end
        if isfield(p.queue, 'remote_recons') && p.queue.remote_recons
            keyboard
        end
        
    end
    
    fprintf('Pausing for 5 seconds.\n')
    fprintf('###############################################################\n\n')
    pause(5);

elseif ~p.getReport.completed
    if isfield(p, 'remote_failed') && p.remote_failed
        try
            io.movefile_fast(fullfile(p.queue.path,'in_progress', p.queue.file_this_recons),fullfile(p.queue.path, 'failed', p.queue.file_this_recons));
        catch
            fprintf('Failed to move file to failed.\n')
        end
    elseif isfield(p.queue, 'file_this_recons')
        try
            fprintf('Reconstruction stopped, moving file to %s.\n', fullfile(p.queue.path, p.queue.file_this_recons))
            io.movefile_fast(fullfile(p.queue.path,'in_progress', p.queue.file_this_recons),fullfile(p.queue.path, p.queue.file_this_recons))
        catch
            fprintf('Failed to move file back to queue search path.\n')
        end
    end
    if isfield(p.queue, 'lockfile') 
        if isempty(p.queue.lockfile)
            if verbose > 2
                p.queue.lockfile = false;
            else
                p.queue.lockfile = true;
            end
        end
        if p.queue.lockfile
            if isempty(p.save_path{1})
                try
                    for ii = 1:length(p.scan_number)
                        p.scan_str{ii} = sprintf(p.scan_string_format, p.scan_number(ii));        % Scan string
                    end
                    p = core.ptycho_prepare_paths(p);
                catch
                    fprintf('Could not find lock file. \n')
                end
            end
            if isfield(p, 'run_name')  && ~isempty(p.run_name)
            for ii=1:length(p.save_path)
                lock_filename = [p.save_path{ii} '/' p.run_name '_lock'];
                if exist(lock_filename, 'file')
                    try
                        delete(lock_filename);
                        fprintf('Removing lock file %s\n',lock_filename)
                    catch
                        fprintf('Removing lock file %s failed\n',lock_filename)
                    end
                end
            end
            end
        end
    end
    if isfield(p, 'remote_file_this_recons')
        try
            fprintf('Removing remote file.\n')
            if exist(p.queue.remote_file_this_recons, 'file')
                delete(p.queue.remote_file_this_recons)
            end
            
            [~, this_file] = fileparts(p.queue.remote_file_this_recons);
            if p.queue.isreplica
                system(['touch ' fullfile(p.queue.remote_path, [this_file '.crash'])]);
            end
            
            this_file = [this_file '.mat'];
            
            if exist(fullfile(p.queue.remote_path, 'in_progress', this_file), 'file')
                delete(fullfile(p.queue.remote_path, 'in_progress', this_file));
            end
            if exist(fullfile(p.queue.remote_path, 'done', this_file), 'file')
                delete(fullfile(p.queue.remote_path, 'done', this_file));
            end
            if exist(fullfile(p.queue.remote_path, 'done', this_file), 'file')
                delete(fullfile(p.queue.remote_path, 'done', this_file));
            end

        catch
            fprintf('Failed to remove remote file.\n')
        end
    end
    
    if ~isempty(p.io.phone_number) && p.io.send_crashed_recon_SMS
        %modified by YJ
        subject = 'Ptycho recon crashed!';
        message = sprintf('User %s''s ptycho recon of %s scan%s crashed on %s', io.get_user_name(), ...
            p.io.data_descriptor, num2str(p.scan_number), io.get_host_name());

        if isfield(p.engines{1},'use_gpu') && p.engines{1}.use_gpu 
            message = strcat(message,'(GPU id ',num2str(p.engines{1}.gpu_id),')');
        end

        io.sendSMS(p.io.phone_number, subject, message);
    end
end

pid = [p.ptycho_matlab_path './utils/.tmp_procID/proc_' num2str(feature('getpid')) '.dat'];
if exist(pid, 'file')
    delete(pid)
end

if ~isempty(p.io.phone_number) && p.io.send_finished_recon_SMS && p.getReport.completed
    %modified by YJ
    subject = 'Ptycho recon completed!';
    message = sprintf('User %s''s ptycho recon of %s scan%s completed on %s', io.get_user_name(), ...
        p.io.data_descriptor, num2str(p.scan_number), io.get_host_name());
    
    if isfield(p.engines{1},'use_gpu') && p.engines{1}.use_gpu 
        message = strcat(message,'(GPU id ',num2str(p.engines{1}.gpu_id),')');
    end
    
    io.sendSMS(p.io.phone_number, subject, message);
end

try
    verbose(struct('prefix', {[]}))
catch
end

end
