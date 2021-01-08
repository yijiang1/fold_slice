%MATLAB_APS - Matlab data preparator for APS hdf5 data
% Writtend by YJ

function p = matlab_aps(p)
import utils.*

%% call prepare data function
for ii=1:length(p.scan_number)
    p.scanID = ii;
    [p] = detector.prep_data.matlab_aps.prepare_data(p);
end

%% combine data structures
data = single([]);
fmask = logical([]);
for ii=1:length(p.scan_number)
    data = cat(3,data,p.detectors(ii).detStorage.data);
    if size(p.detectors(ii).detStorage.fmask,3) > 1
        fmask = cat(3,fmask,logical(round(p.detectors(ii).detStorage.fmask)));
    else
        fmask = cat(3,fmask,repmat(p.detectors(ii).detStorage.fmask,[1 1 size(p.detectors(ii).detStorage.data,3)]));
    end
    p.detectors(ii).detStorage.fmask = [];  % save memory 
    p.detectors(ii).detStorage.data = [];   
end

for ii = 1:p.numobjs
    p.object_size(ii,:) = [size(p.object{ii},1),size(p.object{ii},2)]; 
end

if (isfield(p.detector,'binning')&& p.detector.binning) || (isfield(p.detector,'upsampling')&& p.detector.upsampling)
    p = core.apply_binning(p, 2^(p.detector.binning - p.detector.upsampling));   % modify the p structure after binning 
end

p = detector.prep_data.matlab_ps.prep_data_matlab(p, data, fmask);

% Compare number of points and diffraction patterns
num_difpat = size(p.fmag,3);
verbose(2, 'Number of probe positions: %d', sum(p.numpts));
verbose(2, 'Number of diffraction patterns : %d', num_difpat);
if num_difpat ~= sum(p.numpts)
    error('Number of probe positions (%d) inconsistent with number of diffraction patterns (%d)', sum(p.numpts), num_difpat);
end


end