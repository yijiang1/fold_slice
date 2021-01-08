%EIGER settings for Eiger detector
%   This file is called if you specify p.detector = 'eiger_APS'.
% ** p          p structure
% returns:
% ++ det        detector structure; later accessible via p.detectors(ii).params
%
% see also: detector.load_detector


function [ det ] = medipix3( p )

local_path = fileparts(mfilename('fullpath')); 


%% parameter
det.pixel_size = 55e-6;                 % detector pixel size
det.file_extension = 'h5';              % raw data file extension 
%det.mask = [local_path, '/eiger_valid_mask.mat'];    % detector mask from PSI's detector. Need to change for APS
det.basepath_dir = 'eiger/';            % raw data directory
det.mask_saturated_value = [];          % mask pixels above given value
det.mask_below_value = [];              % mask pixels below given value
det.data_stored = true;                % false == data are generated "onfly", no need to load / store
%det.geometry.sz = [1030 514];           % detector readout size
det.geometry.sz = [1024 1024];           % detector readout size

det.geometry.mask = [];                 % if the mask is larger than the readout 
                                        % geometry (det.geometry.sz), geometry.mask defines the readout for the mask
                                        % e.g. geometry.mask = {[50:500], [50:500]}


det.orientation = [1 0 0];              % [<Transpose> <FlipLR> <FlipUD>]

% additional arguments for image_read
det.image_read_extraargs = {'H5Location', '/eh5/images/'};


%% define directory tree structure
% define a function handle for compiling the scan directory tree. Input
% should be the scan number.
det.read_path = @utils.compile_aps_dirname; 

%% filename
% specify funcion for compiling the filename. Input and output should be p and a 
% temporary structure, containing information about the currently prepared scan. 
% Cf. '+scans/get_filenames_cSAXS.m' 

det.get_filename = @scans.get_filenames_aps;
%{
% additional parameters for the filename function (det.get_filename)
det.filename_pattern{1}.str = '%05d';
det.filename_pattern{1}.del = '_';
det.filename_pattern{1}.content = 'scan';
det.filename_pattern{1}.start = 0;
det.filename_pattern{2}.str = '%05d';
det.filename_pattern{2}.del = '_';
% det.filename_pattern{2}.content = 'burst';
det.filename_pattern{2}.start = 0;
det.filename_pattern{3}.str = '%05d';
det.filename_pattern{3}.del = '.';
% det.filename_pattern{3}.content = 'pos';
det.filename_pattern{3}.start = 0;

if ~isfield(p, 'omny_interferometer') || (strcmp(p.omny_interferometer, 'opos_angle') || strcmp(p.omny_interferometer, 'opos'))
    det.filename_pattern{2}.content = 'burst';
    det.filename_pattern{3}.content = 'pos';
else
    det.filename_pattern{2}.content = 'pos';
    det.filename_pattern{3}.content = 'burst';
end
%}

end

