% function convert(scan, raw_data_path)
% convert .raw files to .hdf5
%  Inputs
%     scan              scan number
%     raw_data_path     path to the raw data (optional)

function convert2hdf5(scan, raw_data_path)

    convertor_path = '/sls/X12SA/data/x12saop/EigerPackage/slsDetectorsPackage/bin/hdf5MakerOMNY'; 
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
           error('RAW data are missing, data cannot be converted')
        else
           delete(sprintf('%s/*.h5',load_dir))
        end
        list_h5 = dir([load_dir, '/run_*.h5']); 
    end
    
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
%            [~,out] = system([cbf_maker_path ' ' fullfile(load_dir,sprintf('run_d0_f00000000%02i000_%i.raw',ind_scans(ii) ,scan))]);     
            systemcall = [convertor_path ' ' fullfile(list_raw(1).folder,list_raw(1).name)];
            fprintf('%s\n',systemcall);
            [~,out] = system(systemcall);
            out
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
            
        h5fileinfo = h5info(fullfile(list_h5.folder,list_h5.name));
            
        nframes_converted = h5fileinfo.Groups.Groups(1).Datasets.Dataspace.Size(3);
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
