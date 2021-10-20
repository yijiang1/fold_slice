%PREPARE_DATA_MATLAB prepares and normalizes data for matlab engines

function [ p ] = prep_data_matlab( p, data, fmask )
import utils.pshift
import utils.verbose
import utils.get_option
import math.fftshift_2D

scanfirstindex = [1 cumsum(p.numpts)+1]; % First index for scan number

% Region of interest [xmin xmax ymin ymax]
% Warning Not implemented for shared scans
% Need to update p.numpts, p.positions, data, fmask, indices
if isfield(p, 'scan') && isfield(p.scan, 'roi') &&  ~isempty(p.scan.roi)
    % Convert to p.positions centered on the object
    if p.share_object
       positions_centered = p.positions - p.object_size(1,:)/2 + p.asize(:)/2;
    else
        for ii = 1:p.numscans
            positions_centered(p.scanidxs{ii},:) = p.positions(p.scanidxs{ii},:) - p.object_size(ii,:)/2;
        end
    end
    xmin = p.scan.roi(1)/p.dx_spec(2);
    xmax = p.scan.roi(2)/p.dx_spec(2);
    ymin = p.scan.roi(3)/p.dx_spec(1);
    ymax = p.scan.roi(4)/p.dx_spec(1);
    % Quick check
    if (xmin>xmax)
        error('ROI is empty, xmax<xmin')
    elseif (ymin>ymax)
        error('ROI is empty, ymax<ymin')
    end
    % Do the comparison and update accordingly [xmin xmax ymin ymax]
    whichtokeep = find(  (positions_centered(:,1) > ymin) & ...
        (positions_centered(:,1)< ymax) & ...
        (positions_centered(:,2) > xmin) & ...
        (positions_centered(:,2) < xmax) );
    % update p.numpts, positions, data, fmask, indices
    data = data(:,:,whichtokeep);
    fmask = fmask(:,:,whichtokeep);
    p.positions = p.positions(whichtokeep,:);
    
    for ii = 1:p.numscans
        %         indaux = find(whichtokeep<scanfirstindex(ii+1),'last'); % Finds the last index of scan ii in the new variables
        indaux = find((whichtokeep>=scanfirstindex(ii))&(whichtokeep<scanfirstindex(ii+1))); % Finds indices of scan ii in the new variables
        if isempty(indaux)
            error('ROI: Scan %s has no points in the ROI. Change the ROI or remove this scan from the scan list',p.scan_str{ii})
        end
        p.numpts(ii) = numel(indaux);
    end
    scanfirstindex = [1 cumsum(p.numpts)+1]; % First index for scan number
    for ii = 1:p.numscans
        p.scanindexrange(ii,:) = [scanfirstindex(ii) scanfirstindex(ii+1)-1];
    end
    verbose(2, 'Region of interest reduced number of points from %d to %d', sum(p.numpts_orig), sum(p.numpts));
    
    % Redo a convenient offset
    p.positions = p.positions - min(p.positions) + 1;
    
   
    % Recompute object sizes
    p = update_object_size(p); 

   
end



% Skip some data point (for testing reduced dose) 
if isfield(p, 'skip_points')
    if ~isempty(p.skip_points) && p.skip_points>1
        offset = mod(p.scan_number, p.skip_points);
        %offset = mod([0 1 2], p.skip_points);
        whichtokeep = [1+offset(1):p.skip_points:p.numpts(1)];
        for idx = 2:length(p.numpts)
            last_idx = sum(p.numpts(1:idx-1));
            whichtokeep = [whichtokeep , (last_idx+1+offset(idx)):p.skip_points:(last_idx+p.numpts(idx))];
        end
        data = data(:,:,whichtokeep);
        fmask = fmask(:,:,whichtokeep);
        p.positions = p.positions(whichtokeep,:);
        for ii = 1:p.numscans
            %         indaux = find(whichtokeep<scanfirstindex(ii+1),'last'); % Finds the last index of scan ii in the new variables
            indaux = find((whichtokeep>=scanfirstindex(ii))&(whichtokeep<scanfirstindex(ii+1))); % Finds indices of scan ii in the new variables
            if isempty(indaux)
                error('ROI: Scan %s has no points in the ROI. Change the ROI or remove this scan from the scan list',p.scan_str{ii})
            end
            p.numpts(ii) = numel(indaux);
        end
        scanfirstindex = [1 cumsum(p.numpts)+1]; % First index for scan number
        for ii = 1:p.numscans
            p.scanindexrange(ii,:) = [scanfirstindex(ii) scanfirstindex(ii+1)-1];
        end
        verbose(2, 'Region of interest reduced number of points from %d to %d', sum(p.numpts_orig), sum(p.numpts));

        % Recompute object sizes
        if p.share_object
            p.object_size = p.asize + max(p.positions,[],1);
            verbose(3, 'Computed object size: %d x %d', p.object_size(1), p.object_size(2));
        else
            for ii = 1:p.numscans
                p.object_size(ii,:) = p.asize + max(p.positions(p.scanindexrange(ii,1):p.scanindexrange(ii,2),:),[],1);
                verbose(3, 'Computed object size: %d x %d', p.object_size(ii,1), p.object_size(ii,2));
            end
        end
    end
end

% % padding
% if any(datasize ~= p.asize)
%     newdata = zeros([p.asize, num_difpat]);
%     offset = floor(.5*(p.asize-datasize));
%     newdata(offset(1) + (1:datasize(1)), offset(2) + (1:datasize(2)), :) = data;
%     data = newdata;
%     clear newdata
%     newfmask = ones([p.asize, size(fmask,3)]);
%     newfmask(offset(1) + (1:datasize(1)), offset(2) + (1:datasize(2)), :) = fmask;
%     fmask = newfmask;
%     clear newfmask
% end


% 'auto_center_data' option is useful if detector shifts provided in
% template cannot be trusted, ie for more than 2 joined scans 
if get_option(p,'auto_center_data') && ~get_option(p,'get_artificial_data')
    for i = 1:p.numscans
        ind = p.scanidxs{i};
        [x,y] = math.center(mean(data(:,:,ind) .* fmask(:,:,ind),3));
        verbose(2,'Auto-shifting diffraction patterns by %i %i px', round(x), round(y))
        data(:,:,ind)  = utils.imshift_fast(data(:,:,ind), x, y, [], 'nearest');
        if size(fmask,3) == size(data,3)
            fmask(:,:,ind) = utils.imshift_fast(fmask(:,:,ind), x, y, [], 'nearest');
        elseif p.numscans == 1
            fmask = utils.imshift_fast(fmask, x, y, [], 'nearest');
        else
           error('Not implemented mask shifting option') 
        end
    end
end
    

p.fmask_per_scan = ndims(fmask) == 3;

% Prepare Fourier projections

% normalization ignores valid mask - should modify
for ii=1:p.numscans
    p.max_sum(ii) = max(sum(sum(data(:,:,p.scanidxs{ii}),1),2),[],3);
end
max_power = max(p.max_sum) / prod(p.asize);
p.renorm = sqrt(1/max_power);

p.Nphot = sum(data(:).*fmask(:)); % Number of photons for regularization normalization (ML optimization)

% store mask pre-fftshifted
p.fmask = fftshift_2D(fmask);
clear fmask 
% precalculate modulus of data, normalize and fftshift
p.fmag = fftshift_2D(sqrt(data)) * p.renorm;


end

