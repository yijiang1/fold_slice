%CONVERT2HDF5_WRAPPER converts Eiger 1.5M raw data files to HDF5 and
% deletes the raw files if the conversion has finished successfully
% convert2hdf5_wrapper(raw_data_path)
% 
% ** raw_data_path      path to the eiger directory, e.g. ~/Data10/
%
% *optional*
% ** scanID             start at the given scan number
% 
% EXAMPLES:
%   % start at scan number 1:
%       convert2hdf5_wrapper('~/Data10/');
%
%   % start at scan number 150:
%       convert2hdf5_wrapper('~/Data10/', 150);     
%
% Pleas note that the script is designed to be used during an ongoing
% measurement, and therefore only converts n-1 datasets, that is it waits
% until the next measurement has started.

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
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


function convert2hdf5_wrapper(raw_data_path, varargin)
import utils.*
if nargin > 1
    scanID = varargin{1};
else
    scanID = 1;
end
    while true
        [started, newScan, specDatFile] = beamline.next_scan_started(raw_data_path, scanID);
        if started
            convert2hdf5(scanID, raw_data_path, specDatFile);
            fprintf('Converting scan %d\n', scanID);
            scanID = newScan;
        else
            fprintf('Waiting for next scan to start.\n');
            pause(1);
        end
    end
        

end


function convert2hdf5(scan, raw_data_path, specDatFile)

    % some defaults
    convertor_path = '~/Data10/bin/eiger1p5M_converter/hdf5MakerOMNY'; 
    xmlLayoutFile = '~/Data10/bin/nexus/layout.xml';
    orchestraPath = '~/Data10/specES1/scan_positions/';
    specParser = '~/Data10/matlab/+io/spec_reader/spec_reader';
    
    % check if orchestraPath exists
    if exist(orchestraPath, 'dir')
        orchestraPath = ['--orchestra '  orchestraPath];
    else
        orchestraPath = '';
    end
    
    load_dir = utils.compile_x12sa_dirname(scan);
    if exist('raw_data_path','var')&&exist(fullfile(raw_data_path,load_dir),'dir')
        load_dir = fullfile(raw_data_path,load_dir);
    elseif exist(['~/Data10/eiger_4/'],'dir')
        load_dir = ['~/Data10/eiger_4/' load_dir]; 
    elseif exist([raw_data_path,'eiger_4/'])
        load_dir = [raw_data_path,'/eiger_4/' load_dir]; 
    elseif exist([raw_data_path,'/eigeromny/'])
        load_dir = [raw_data_path,'/eigeromny/' load_dir]; 
    end
    
    if ~exist(load_dir, 'dir')
        warning('Raw data path %s not found', load_dir)
        return
    end

    testDir = [load_dir, '/deleteMe'];
    % test for write permissions by creating a folder and then deleting it
    isWritable = mkdir(testDir);
    % check if directory creation was successful
    if isWritable == 1
        rmdir(fullfile(testDir));
    end

    list_h5 = dir([load_dir, '/run_*.h5']); 
        
    
    file_sizes = [list_h5.bytes]; 
    if any(file_sizes < 1e6)  % find files < 1MB
        warning('H5 files in scan %i seem damaged, generate again ... ', scan)
        list_raw = dir([load_dir, '/run_d0_f0000000*.raw']);
        if isempty(list_raw)
           warning('RAW data is missing, data cannot be converted')
           return
        else
           delete(sprintf('%s/*.h5',load_dir))
        end
        list_h5 = dir([load_dir, '/run_*.h5']); 
    end
%     toc
    if isempty(list_h5)    
        if ~isWritable
            warning('Conversion failed because folder %s is not writable', load_dir)
            return
        end
        
        list_raw = dir(fullfile(load_dir, 'run_d0_f0000000*.raw')); 

        Nscans = length(list_raw); 
        
        for ii = 1:Nscans
            ind_scans(ii) = str2num(list_raw(ii).name(16:17)); 
        end
        
        for ii = 1:Nscans
            systemcall = [convertor_path ' ' fullfile(list_raw(1).folder,list_raw(1).name)];
            fprintf('%s\n',systemcall);
            [stat,out] = system(systemcall);
            systemcall = sprintf('%s -s %s --scanNr %u --hdf5 --xmlLayout %s -o %s %s', specParser, specDatFile, scan, xmlLayoutFile, fullfile(load_dir, sprintf('run_%05d_000000000000.h5',scan)), orchestraPath);
            [stat, out_spec] = system(systemcall);
        end
        
        list_h5 = dir([load_dir, '/*.h5']);
        if isempty(list_h5)
           error(sprintf('After conversion did not find any h5 in %s\n',load_dir))
           return
        end
        if numel(list_h5)>1
           error(sprintf('After conversion I found more than one h5 in %s\n',load_dir))
           return
        end
            
        h5fileinfo = h5info(fullfile(list_h5.folder,list_h5.name), '/entry/instrument/eiger_4/data');
            
        nframes_converted = h5fileinfo.Dataspace.Size(3);
        out = splitlines(out); 
        nframes_expected = str2num(out{end-2}(14:end)); 
        fprintf('Frames expected: %i, frames converted %i \n', nframes_expected, nframes_converted)
        if nframes_converted == nframes_expected
            fprintf('Scan %i succefully converted to H5\n', scan); 
            delete(sprintf('%s/*.raw',load_dir))
        else
            error('Scan %i WAS NOT CONVERTED to H5\n', scan)    
            delete(sprintf('%s/*.h5',load_dir))
        end
    end
end
