%APPLY_CIRCULAR_MASK applies a circular mask to fmask, which will be
%applied to data later
%
% ** p      p structure
% returns:
% ++ p      updated p structure
%
% see also: detector.prep_data.matlab_ps.matlab_aps


function [ p ] = apply_circular_mask( p )

if isfield(p.detector,'circ_mask') && p.detector.circ_mask>0
    utils.verbose(2, 'Apply a circular mask (radius=%d pixels) to diffraction patterns',p.detector.circ_mask)
    fmask = p.detectors(p.scanID).detStorage.fmask;
    circ_mask = utils.make_circular_mask([size(fmask,1),size(fmask,1)], p.detector.circ_mask);
    p.detectors(p.scanID).detStorage.fmask = fmask .*circ_mask;
end

end

