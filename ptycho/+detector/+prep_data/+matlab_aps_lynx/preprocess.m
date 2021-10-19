%PREPROCESS preprocess functions
%   e.g. select area before reading from disk
%
% ** p      p structure
%
% returns:
% ++ p      p structure
%
% see also: detector.prep_data.matlab_ps.prepare_data
%

function [ p ] = preprocess( p )

detStorage = p.detectors(p.scanID).detStorage;


p.detectors(p.scanID).params.image_read_extraargs{end+1} = 'RowFrom';
p.detectors(p.scanID).params.image_read_extraargs{end+1} = detStorage.lim_inf(1);
p.detectors(p.scanID).params.image_read_extraargs{end+1} = 'RowTo';
p.detectors(p.scanID).params.image_read_extraargs{end+1} = detStorage.lim_sup(1);
p.detectors(p.scanID).params.image_read_extraargs{end+1} = 'ColumnFrom';
p.detectors(p.scanID).params.image_read_extraargs{end+1} = detStorage.lim_inf(2);
p.detectors(p.scanID).params.image_read_extraargs{end+1} = 'ColumnTo';
p.detectors(p.scanID).params.image_read_extraargs{end+1} = detStorage.lim_sup(2);

if ~p.prealign_FP
    detStorage.read_size = p.asize;
else
    detStorage.read_size = p.prealign.asize;
end


end

