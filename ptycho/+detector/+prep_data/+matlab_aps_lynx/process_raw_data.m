%PROCESS_RAW_DATA postprocess functions
%   applied to the 4D raw data
%
% ** p      p structure
% returns:
% ++ p      updated p structure
%
% see also: detector.prep_data.matlab_ps.matlab_ps


function [ p ] = process_raw_data( p )

detStorage = p.detectors(p.scanID).detStorage;

if p.prealign_FP
    p = core.FPM.FP_prealign(p);
    detStorage.lim_inf = detStorage.ctr-p.asize/2;
    detStorage.lim_sup = detStorage.ctr+p.asize/2-1;
end

% sum up data from burst scans
if p.scan.is_cont
    detStorage.data = sum(detStorage.data,4);
else
    sz = size(detStorage.data);
    detStorage.data = reshape(detStorage.data, sz(1), sz(2), []);
end


end

