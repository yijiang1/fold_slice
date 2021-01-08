%PREPARE_DATA prepare the raw data for the ptychographic reconstruction
%   
% ** p          p structure
%
% returns:
% ++ p          p structure
%
%   *Note:*
%   All (default) functions are located in +detector/+prep_data/+matlab_ps. If You want
%   to add a new detector, please create a new package directory
%   (+detector_name) with a parameter file detector_name.m. Functions in a
%   package detectory directory will overload similar functions in the default directory.
%
% see also: detector.prep_data.matlab_ps.matlab_ps
%
%


function [ p ] = prepare_data( p )
import utils.*
import io.image_read

verbose(1, 'Preparing data using matlab APS data preparation.')


%% initialize detector
scanID = p.scanID;
detStorage = p.detectors(scanID).detStorage;

if size(p.ctr,1) == 1
    center = p.ctr; % save centers for all diffraction patterns
else
    center = p.ctr(scanID,:); % different centers for different diffraction patterns
end

check_ctr = 'auto';
if ~isempty(center)
    check_ctr = 'inter';
    if ~strcmp(center, 'inter')
        detStorage.ctr = center;
        verbose(2,['Using supplied center: ', num2str(detStorage.ctr)]);
        check_ctr = 'no';
    end
end
detStorage.check_ctr = check_ctr; 


%% load mask
if isfield(p.detectors(scanID).params, 'mask') && ~isempty(p.detectors(scanID).params.mask)
    load(p.detectors(scanID).params.mask);
    detStorage.mask = logical(mask);
elseif p.detectors(scanID).params.data_stored 
    detStorage.mask = ones(p.detectors(scanID).params.geometry.sz, 'logical');
else
    detStorage.mask = logical([]); 
end


%% crop mask to readout geometry (if necessary)

[p] = p.detectors(scanID).funcs.mask_geometry(p);

%% get center and readout size

[p] = p.detectors(scanID).funcs.get_center(p);

%% preprocess

[p] = p.detectors(scanID).funcs.preprocess(p);

%% load data

[p] = p.detectors(scanID).funcs.load_data(p);

%% process raw data

[p] = p.detectors(scanID).funcs.process_raw_data(p);

%% mask saturated values

[p] = p.detectors(scanID).funcs.get_mask(p);

%% apply a circular mask to diffraction patterns. Added by YJ.

[p] = p.detectors(scanID).funcs.apply_circular_mask(p);

%% apply binning on the measured data / mask 
[p] = p.detectors(scanID).funcs.binning(p);

%% final step postprocessing, e.g. background subtraction 

[p] = p.detectors(scanID).funcs.postprocess(p);


end

