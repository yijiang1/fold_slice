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
    if ~par.use_OMNY_file_angles
        S = io.spec_read(par.base_path,'ScanNr',scans);
        for ii = 1:Nscans
            angles(ii)=S{ii}.samroy;
        end
    else
        [S, errflag] = beamline.read_omny_angles(par.OMNY_angle_file,scans, tomo_id);
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
