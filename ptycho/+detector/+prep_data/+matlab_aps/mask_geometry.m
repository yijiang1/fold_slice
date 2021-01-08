%MASK_GEOMETRY crop mask to fit defined readout geometries
%   if no geometry is specified, the full mask is used
% ** p      p structure
%
% returns:
% ++ p      p structure
%
% see also: detector.prep_data.matlab_ps.prepare_data
%

function [ p ] = mask_geometry( p )

scanID = p.scanID;
detStorage = p.detectors(scanID).detStorage;

if isfield(p.detectors(scanID).params, 'geometry') && ~isempty(p.detectors(scanID).params.geometry.mask)
    detStorage.mask = detStorage.mask(p.detectors(ii).params.geometry.mask{1},p.detectors(ii).params.geometries.mask{2},:);
end

end

