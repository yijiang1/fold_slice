%EIGER1P5M settings for Eiger 1.5M detector
%   This file is called if you specify p.detector = 'eiger1p5M'
%
% ** p          p structure
% returns:
% ++ det        detector structure; later accessible via p.detectors(ii).params
%
% see also: detector.load_detector

function [ det ] = eiger1p5M_APS( p )

local_path = fileparts(mfilename('fullpath')); 


%% parameter
det.pixel_size = 75e-6;                 % detector pixel size
det.file_extension = 'h5';             % raw data file extension 
det.mask = [local_path, '/eiger1p5m_APS_mask.mat'];    % detector mask
if exist(fullfile(p.raw_data_path{1}, 'eiger_4'))
    det.basepath_dir = 'eiger_4/';           % raw data directory
else
    det.basepath_dir = 'eigeromny/';           % raw data directory, legacy, remove in future 
end
det.mask_saturated_value = [];          % mask pixels above given value
det.mask_below_value = [];              % mask pixels below given value
det.data_stored = true;               % false == data are generated "onfly", no need to load / store
det.detposmotor = 'dettr';              % name of spec motor used for checking detector position
det.geometry.sz = [1030 1614];          % detector readout size
det.geometry.mask = [];                 % if the mask is larger than the readout 
                                        % geometry (det.geometry.sz), geometry.mask defines the readout for the mask
                                        % e.g. geometry.mask = {[50:500], [50:500]}


% det.orientation = [0 0 1];              % [<Transpose> <FlipLR> <FlipUD>]
% changed from 30/08/2018
det.orientation = [1 0 0];              % [<Transpose> <FlipLR> <FlipUD>]


% additional arguments for image_read
if strcmpi(det.file_extension, 'h5')
    det.image_read_extraargs = {'H5Location', '/entry/instrument/eiger_4/data/'};
else
    det.image_read_extraargs = {};
end

%% define directory tree structure
% define a function handle for compiling the scan directory tree. Input
% should be the scan number.
det.read_path = @utils.compile_x12sa_dirname; 

%% filename
% specify funcion for compiling the filename. Input and output should be p and a 
% temporary structure, containing information about the currently prepared scan. 
% Cf. '+scans/get_filenames_cSAXS.m' 

%det.get_filename = @scans.get_filenames_cSAXS;
det.get_filename = @scans.get_filenames_aps_lynx;

% additional parameters for the filename function (det.get_filename)
% define filename pattern
% det.read_path_format = '%05d/';
% det.filename_pattern{1}.str = '%05d';
% det.filename_pattern{1}.del = '_';
% det.filename_pattern{1}.content = 'scan';
% det.filename_pattern{1}.start = 1;
% det.filename_pattern{2}.str = '%05d';
% det.filename_pattern{2}.del = '.';
% det.filename_pattern{2}.content = 'pos';
% det.filename_pattern{2}.start = 1;
% 
% det.data_prefix = 'run_';


% for scan = p.scan_number
% %     if ~isempty(dir(fullfile(p.raw_data_path{1}, 'eiger_4', utils.compile_x12sa_dirname(scan), '*.cbf')))
% %         det.file_extension = 'cbf';
% %     end
%     %%%% convert data from raw to cbf or hdf5
%     switch det.file_extension
%         case 'cbf'
%             detector.eiger1p5M.convert(scan, p.raw_data_path{1});
%         case 'h5'
%             detector.eiger1p5M.convert2hdf5(scan, p.raw_data_path{1});
%         otherwise
%             error('Unknown file extension for Eiger 1.5M.')
%     end
% 
% end


end

