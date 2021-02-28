%   LOAD_ANGLES_APS  load tomopgrahy angles for given scan numbers 
%   Directly load from original master .h5 files
%   created by YJ based on PSI's function

%   [par,  angles] = load_angles_aps(par, scans, plot_angles)
%   Inputs: 
%       **par               tomo parameter structure 
%       **scans           - list of loaded scan numbers 
%       **tomo_id         - indetification number of the sample, default = []
%       **plot_angles     - plot loaded angles, default == true 
%   *returns*
%       ++par               tomo parameter structure 
%       ++angles            loaded angles 

function [par, angles] = load_angles_aps(par, scans, plot_angles)
    warning on
    Nscans = length(scans);
    angles = nan(Nscans,1);
    hasAngle = ones(Nscans,1);
    if isfield(par,par.angle.filesuffix) && ~isempty(par.angle.filesuffix)
        file_suffix = par.angle.filesuffix;
    else
        file_suffix = '_master.h5'; %default for velociprobe data outputs
    end
    
    if isfield(par.angle,'h5path') && ~isempty(par.angle.h5path)
        h5path = par.angle.h5path;
    else
        h5path = '/entry/sample/goniometer/chi_start';
    end
    
    %%
    wb = waitbar(0,'1','Name','Loading ptycho-tomo projection angles...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(wb,'canceling',0);
    for i=1:Nscans
         % Check for clicked Cancel button
        if getappdata(wb,'canceling')
            break
        end

        filename = strcat(par.base_path, 'ptycho/',sprintf(par.scan_string_format, scans(i)),'/',sprintf(par.scan_string_format, scans(i)),file_suffix);
        if ~isempty(filename)
            try
            	angle_temp = h5read(filename,h5path);
                angles(i) = angle_temp(1);
                status = [sprintf(par.scan_string_format, scans(i)), ' angle = ',num2str(angles(i))];
            catch
                disp(['Reading angle failed for ', sprintf(par.scan_string_format, scans(i))]);
                disp(strcat('Check angle h5path:',h5path))
                disp(filename)
                status = ['Reading angle failed for ', sprintf(par.scan_string_format, scans(i))];
            end
        else
            hasAngle(i) = 0;
            disp(['No angle found for ',sprintf(par.scan_string_format, scans(i))])
            status = ['No angle found for ',sprintf(par.scan_string_format, scans(i))];
        end
        
        % Update waitbar and message
        %waitbar(i/Nscans,wb,sprintf(par.scan_string_format, scans(i)))
        waitbar(i/Nscans,wb,status)
        
    end
    delete(wb)

    
    % legacy code - read angles from processed h5 files
    %{ 
    for i=1:Nscans
        file = find_projection_files_names_aps(par, scans(i));
        if ~isempty(file)
            angles(i) = h5read(file,'/angle');
        else
            hasAngle(i) = 0;
            disp(strcat('No angle found for scan ',num2str(scans(i))))
        end
    end
    %}
    %% process angles
    % remove scan without angle
    angles = angles(hasAngle==1);
    scans = scans(hasAngle==1);

    % remove duplicted scan numbers 
    [~,ind] = unique(scans, 'last');   % take the !last! occurence of the scan, assume that the second measurement was better 
    angles = angles(ind);    % Angles not repeated in scan
    scans = scans(ind);
    %subtomos = subtomos(ind);

    % take only unique angles, measure uniqueness
    if par.remove_duplicated_angles 
        [~,ind] = unique(angles, 'last');   % take the !last! occurence of the angle, assume that the second measurement was better 
        if length(angles) ~= length(ind)
            warning('Removed %i duplicated angles', length(angles) - length(ind))
        end
    else
        [~,ind] = sort(angles);     
    end
    angles = angles(ind);    % Angles not repeated in scan
    scans = scans(ind);
    %subtomos = subtomos(ind);

    angles = angles + 0.1;  % avoid the angles to be too well aligned with pixels, ie avoid exact angles 0, 90, 180, ... 

    par.scanstomo = scans;
    %par.subtomos = subtomos;

    par.num_proj=numel(par.scanstomo);
    [anglessort,indsortangle] = sort(angles);
    
    if par.sort_by_angle
        angles = angles(indsortangle);
        par.scanstomo = par.scanstomo(indsortangle);
        %par.subtomos = par.subtomos(indsortangle);
    else  % sort by scan number 
        [~,indsortscan] = sort( par.scanstomo);
        angles = angles(indsortscan);
        par.scanstomo = par.scanstomo(indsortscan);
        %par.subtomos = par.subtomos(indsortscan);
    end

    if par.verbose_level && plot_angles
        plotting.smart_figure(1);
        subplot(2,1,1)
        plot(par.scanstomo,angles,'ob'); grid on; 
        xlim(par.scanstomo([1,end]))
        %legend('Tilt angles')
        xlabel('Scan #')
        ylabel('Tilt angles')

        subplot(2,1,2)
        plot(diff(anglessort))
        ylabel('Angle increment')
        %title('Angular spacing'); 
        grid on; 
        xlim([1,par.num_proj-1])
        if par.windowautopos
            screensize = get( groot, 'Screensize' );
            win_size = [946 815]; 
            set(gcf,'Outerposition',[139 min(163,screensize(4)-win_size(2))  win_size]);  %[left, bottom, width, height]
        end
        title('Measured angles')
        drawnow 
    end
end
