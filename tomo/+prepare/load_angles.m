%   LOAD_ANGLES  load tomopgrahy angles for given scan numbers or tomo_id 
%
%   [par,  angles] = load_angles(par, scans, tomo_id, plot_angles)
%   Inputs: 
%       **par               tomo parameter structure 
%       **scans           - list of loaded scan numbers 
%       **tomo_id         - indetification number of the sample, default = []
%       **plot_angles     - plot loaded angles, default == true 
%   *returns*
%       ++par               tomo parameter structure 
%       ++angles            loaded angles 

function [par,  angles] = load_angles(par, scans, tomo_id, plot_angles)
    if nargin < 4
        plot_angles = true; 
    end
    warning on
    if nargin < 3
        tomo_id = [];
    end
    Nscans = length(scans);
    angles = nan(Nscans,1);
    if isfield(par,'scannum_offset')
        scannum_offset = par.scannum_offset;
    else
        scannum_offset = 0;
    end
    if ~par.use_OMNY_file_angles
        S = io.spec_read(par.base_path,'ScanNr',scans);
        for ii = 1:Nscans
            angles(ii)=S{ii}.samroy;
        end
    else
        [S, errflag] = beamline.read_omny_angles(par.OMNY_angle_file,scans, tomo_id, scannum_offset);
        if errflag
            disp(['Not all scans found in ' par.OMNY_angle_file])
            disp(['I will remove the angles not found and show you some plots anyway'])
        end
        angles=S.readout_angle(:).';
        scans = S.scan(:).';
        subtomos = S.subtomo_num(:).';
        if isfield(S,'tomo_id')
            if any(S.tomo_id ~= S.tomo_id(1))
                warning('tomo_id number is not the same for all scans')
            end
            par.tomo_id = unique(S.tomo_id);
        else
            par.tomo_id = [] ; 
        end
        par.sample_name = S.sample_name{1};
    end
   
    % remove duplicted scan numbers 
    [~,ind] = unique(scans, 'last');   % take the !last! occurence of the scan, assume that the second measurement was better 
    angles = angles(ind);    % Angles not repeated in scan
    scans = scans(ind);
    subtomos = subtomos(ind);


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
    subtomos = subtomos(ind);

    if isfield(par,'angle_offset') && par.angle_offset ~=0 
        angles = angles + par.angle_offset;  % avoid the angles to be too well aligned with pixels, ie avoid exact angles 0, 90, 180, ... 
    end
    
    par.scanstomo = scans;
    par.subtomos = subtomos;

    par.num_proj=numel(par.scanstomo);
    [anglessort,indsortangle] = sort(angles);
    
    if par.sort_by_angle
        angles = angles(indsortangle);
        par.scanstomo = par.scanstomo(indsortangle);
        par.subtomos = par.subtomos(indsortangle);
    else  % sort by scan number 
        [~,indsortscan] = sort( par.scanstomo);
        angles = angles(indsortscan);
        par.scanstomo = par.scanstomo(indsortscan);
        par.subtomos = par.subtomos(indsortscan);
    end

    if par.verbose_level && plot_angles
        plotting.smart_figure(1);
        subplot(2,1,1)
        plot(par.scanstomo,angles,'ob'); grid on; 
        %par.scanstomo(1)
        %par.scanstomo(end)
        xlim(par.scanstomo([1,end]))
        legend('Spec angles')
        xlabel('Scan #')
        subplot(2,1,2)
        plot(diff(anglessort))
        title('Angular spacing'); grid on; 
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
