%LOAD_PREPARED_DATA Load prepared data file and convert it into the default Matlab
%structure
%
%   filename             path and filename of the h5 file
%
%   *optional*
%   return_intensity     return intensity or magnitude; default false (= return magnitude)
%   scan                 select scan, either integer or array
%   enum                 return only selected frames 
%   return_fftshifted    return results fftshifted, default == true 
%
%   *returns*: 
%   fmag                fourier magnitudes of the measured data, ie fftshift(sqrt(data))
%   fmask               mask of the fourier magnitudes, 1 for bad pixels, 0 for other
%   pos                 scanning positions (Npos x 2 array)
%   max_power           maximal intesity (max(sum(sum(fmag,1),2),[],3) / numel(fmag(:,:,1));)
%   scanindexrange      indices corresponding to each of the scans 
%   max_sum             something stored in h5_data.measurement.(['n' num2str(ii-1)]).Attributes.max_sum;
%
%   Examples:
%       [fmag, fmask, ~] = load_prepared_data('~/Data10/analysis/S00668/S00668_S00669_data_400x400.h5');
%       [fmag, fmask, pos] = load_prepared_data('~/Data10/analysis/S00668/S00668_S00669_data_400x400.h5');
%
%       % load intensities
%       [I, ~, ~] = load_prepared_data('~/Data10/analysis/S00668/S00668_S00669_data_400x400.h5', true);
%   
%       % load data from second scan
%       [fmag, fmask, pos] = load_prepared_data('~/Data10/analysis/S00668/S00668_S00669_data_400x400.h5', false, 2);
%
%       % load data from scan 1 and 3
%       [fmag, fmask, pos] = load_prepared_data('~/Data10/analysis/S00668/S00668_S00669_data_400x400.h5', false, [1 3]);


%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
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

function [ fmag, fmask, pos, max_power, scanindexrange, max_sum ] = load_prepared_data( filename, return_intensity, enum, return_fftshifted )
import io.HDF.hdf5_load

if nargin < 2
    return_intensity = false;  % return intensity as measured by detector 
end
if nargin < 3
    enum = [];                 % return only selected frames 
end
if nargin < 4
    return_fftshifted = true;  % return data fftshifted, !! DEFAULT == true !! 
end

    
if ~exist(filename, 'file'); error('Cannot load %s', filename); end

%% compatility to load also old matlab datasets 
[~,~,ext]=fileparts(filename); 
if strcmpi(ext, '.mat')
    d = load(filename); 
    fmag = d.data; 
    fmask = d.fmask; 
    if ~return_intensity
        % normalize to provide similar data as in h5-mex
        max_power = max(sum(sum(fmag,1),2),[],3) / numel(fmag(:,:,1));
        renorm = sqrt(1/max_power);
        fmag = sqrt(fmag)*renorm;
    end
    if return_fftshifted
        fmag = math.fftshift_2D(fmag); 
        fmask = math.fftshift_2D(fmask); 
    end
    return
end
    
    
%% load h5 file
recon  = false;
inf = h5info(filename);
for ii=1:numel(inf.Groups)
    if strcmpi(inf.Groups(ii).Name, '/reconstruction')
        recon = true;
    end
end

if recon
    h5_data = hdf5_load(filename, '/measurement/data', '-a');
else
    h5_data = hdf5_load(filename, '-a');
end

if isfield(h5_data, 'reconstruction')
    h5_data = h5_data.measurement.data;
end

%% check hdf5 data verion (python or mex data prep?)
if isfield(h5_data, 'measurements')
    h5_version = 'mex';
else
    h5_version = 'LibDetXR';
end

switch h5_version
    case 'mex'
        warning('Outdated data format.')
        asize = size(h5_data.measurements.measurement_0.diff_pat.Value);
        fn = fieldnames(h5_data.measurements);
        numpts = length(fn)-1;
        fmag = zeros(asize(1), asize(2), numpts);
        fmask = ones(asize(1), asize(2), numpts);
        pos = zeros(numpts, 2);
        
        if isempty(enum) 
            enum = 1:length(fieldnames(h5_data.detectors))-1;
        end
        
        for ii=enum 
            
            fmaskdet{ii} = zeros(asize(1),asize(2));
            modules = transpose(h5_data.detectors.(['detector_' num2str(ii-1)]).modules.Value);
            
            numrows = modules(:,1);
            numcols = modules(:,2);
            indbeginmody = modules(:,3)+1;
            indbeginmodx = modules(:,4)+1;
            indendmody = numrows + indbeginmody -1;
            indendmodx = numcols + indbeginmodx -1;
            nummody = length(indbeginmody);
            nummodx = length(indbeginmodx);
            
            for kk = 1:nummody
                for jj = 1:nummodx
                    fmaskdet{ii}(indbeginmody(kk):indendmody(kk),indbeginmodx(jj):indendmodx(jj))=1;
                end
            end
        end
        scanindexrange = ones([2 length(enum)]);
        cid=1;
        for ii=1:numpts
            cf =['measurement_' num2str(ii-1)];
            det = h5_data.measurements.(cf).Attributes.detector;
            if any(enum==det+1)
                if ii>scanindexrange(2,det+1)
                    scanindexrange(2,det+1) = ii;
                end

                fmag(:,:,cid) = transpose(h5_data.measurements.(cf).diff_pat.Value);
                if isfield(h5_data.measurements.(cf), 'bad_pixels')
                    bp_temp = h5_data.measurements.(cf).bad_pixels.Value;
                else
                    bp_temp = [];
                end
                pos(cid, :) = h5_data.measurements.(cf).Attributes.position;
                fmask(:,:,cid) = fmaskdet{det+1};
                if ~isempty(bp_temp)
                    for kk=1:size(bp_temp,2)
                        fmask(bp_temp(1,kk)+1,bp_temp(2,kk)+1,cid) = 0;
                    end
                end
                cid = cid + 1;
            end
            
        end
        fmag = fmag(:,:,1:cid-1);
        fmask = fmask(:,:,1:cid-1);
        pos = pos(1:cid-1,:);
%         attr = hdf5_load(filename, '/measurements/','-sa');
        max_power = h5_data.measurements.Attributes.max_power;
        fmask = logical(fmask);
        if return_intensity
            renorm = sqrt(1/max_power);
            fmag = (fmag/renorm).^2;
        end
        
        scanindexrange = scanindexrange';
        for ii=2:length(enum)
            scanindexrange(ii,1) = scanindexrange(ii-1,2)+1;
        end



            
        
        
    case 'LibDetXR'
        
        %% get data dims
        fmag_dim(1) = 0;
        fmag_dim(2) = 1;
        if isempty(enum) 
            enum = 1:length(fieldnames(h5_data.measurement))-1;
        end
        
        for ii=enum
            fmag_temp{ii} = h5_data.measurement.(['n' num2str(ii-1)]).data.Value;
            fmag_dim(ii+2) = size(fmag_temp{ii},3);
        end
        
        asize = size(fmag_temp{enum(1)});
        if asize(1) ~= asize(2)
            error('Loading of asymmetric prepated datasets not supported, use p.force_preparation_data=true')
        end
        
        fmag = zeros(asize(1), asize(2), sum(fmag_dim)-1, 'single');
        fmask = ones(asize(1), asize(2), sum(fmag_dim)-1, 'logical');
        pos = zeros(sum(fmag_dim)-1, 2);
        max_sum = zeros(length(enum), 1);
        
        %% load modules to prepare fmask and load everything into containers
        if isempty(enum) 
            enum = 1:length(fieldnames(h5_data.detector));
        end
        
        for ii=enum        
            % get modules for mask
            fmaskdet{ii} = zeros(asize(1),asize(2), 'logical');
            modules = transpose(h5_data.detector.(['n' num2str(ii-1)]).modules.Value);
            
            numrows = modules(:,1);
            numcols = modules(:,2);
            indbeginmody = modules(:,3)+1;
            indbeginmodx = modules(:,4)+1;
            indendmody = numrows + indbeginmody -1;
            indendmodx = numcols + indbeginmodx -1;
            nummody = length(indbeginmody);
            nummodx = length(indbeginmodx);
            
            for kk = 1:nummody
                for jj = 1:nummodx
                    fmaskdet{ii}(indbeginmody(kk):indendmody(kk),indbeginmodx(jj):indendmodx(jj))=1;
                end
            end
            temp_range = sum(fmag_dim(1:ii+1)):sum(fmag_dim(1:ii+2))-1;
            fmask(:,:,temp_range) = repmat(fmaskdet{ii},[1,1,fmag_dim(ii+2)]);
            
            % get bad pixels
            if isfield(h5_data.detector.(['n' num2str(ii-1)]), 'bad_pixels')
                bp = h5_data.detector.(['n' num2str(ii-1)]).bad_pixels.Value;
                for kk=1:size(bp,2)
                    fmask(bp(1,kk)+1,bp(2,kk)+1,temp_range) = 0;
                end
            else
                if isfield(h5_data.measurement.(['n' num2str(ii-1)]), 'bad_pixels')
                    bp = h5_data.measurement.(['n' num2str(ii-1)]).bad_pixels.Value;
                    bpi = h5_data.measurement.(['n' num2str(ii-1)]).bad_pixels_index.Value;
                    assert(length(bpi)==length(temp_range), 'Number of frames does not match the number of bad pixel datasets.')
                    offset = 1;
                    for kk=1:length(bpi)
                        for jj=offset:bpi(kk)
                            fmask(bp(1,jj)+1, bp(2,jj)+1, temp_range(kk)) = 0;
                        end
                        offset = bpi(kk);
                    end
                end
            end
            fmag(:,:,temp_range) = permute(fmag_temp{ii}, [2 1 3]);
            pos_temp = h5_data.measurement.(['n' num2str(ii-1)]).positions.Value;
            pos(temp_range,1) = pos_temp(1,:);
            pos(temp_range,2) = pos_temp(2,:);
            
            max_sum(enum) = h5_data.measurement.(['n' num2str(ii-1)]).Attributes.max_sum;

            
        end
        
        
        % normalize to provide similar data as in h5-mex
        max_power = max(sum(sum(fmag,1),2),[],3) / numel(fmag(:,:,1));
        renorm = sqrt(1/max_power);
        if ~return_intensity
            fmag = sqrt(fmag)*renorm;
        end
        fmask = logical(fmask);
        
        scanindexrange = zeros(numel(enum),2);
        scanindexrange(1,:) = fmag_dim(2:3);
        for ii=2:numel(enum)
            scanindexrange(ii,1) = scanindexrange(ii-1,2)+1;
            scanindexrange(ii,2) = scanindexrange(ii-1,2)+fmag_dim(ii+2);
        end
        
    otherwise
        error('Unknown HDF5 data structure!')
end

if ~return_fftshifted
    % return data as seen by detector, 
    fmag = math.ifftshift_2D(fmag); 
    fmask = math.ifftshift_2D(fmask);
end


end

