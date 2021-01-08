%GET_FILENAMES_APS compile filenames of APS' hdf5 data files
% written by YJ based on PSI's code

function [p] = get_filenames_aps(p)
import utils.*

% get detector paramters
det = p.detectors(p.scanID).params;
read_path = p.raw_data_path_full{p.scanID};
detStorage = p.detectors(p.scanID).detStorage;

detStorage.files = [];

for ii = 1:p.numscans
    detStorage.files{ii} = strcat(p.base_path,sprintf(p.scan.format, p.scan_number(ii)),'/data_roi',p.scan.roi_label,'_dp.hdf5');   
end

%{
if isfield(det, 'filename_pattern')
    if iscell(det.filename_pattern)
        %% if filename patterns exist, use them to restrict the file search
        det.filename_pattern_full = [p.detector.data_prefix];
        fill = [];
        
        for ii=1:size(det.filename_pattern,2)
            fill = [fill det.filename_pattern{ii}.str det.filename_pattern{ii}.del];
            det.filename_pattern_full = [det.filename_pattern_full det.filename_pattern{ii}.str det.filename_pattern{ii}.del];
            del{ii} = det.filename_pattern{ii}.del;
            
            switch det.filename_pattern{ii}.content
                case 'burst'
                    burst = ii;
                case 'scan'
                    scan = ii;
                case 'pos'
                    pos = ii;
            end
            
        end
        
        det.filename_pattern_full = [det.filename_pattern_full det.file_extension];
        
        if exist('burst', 'var')
            det.filename_pattern_burst = [p.detector.data_prefix];
            parts = strsplit(fill, del);
            parts{burst} = '*';
            det.filename_pattern_burst = [det.filename_pattern_burst strjoin(parts, del) det.file_extension];
        end
        
        if exist('pos', 'var')
            det.filename_pattern_pos = [p.detector.data_prefix];
            parts = strsplit(fill, del);
            parts{pos} = '*';
            det.filename_pattern_pos = [det.filename_pattern_pos strjoin(parts, del) det.file_extension];
        end
        
        if exist('scan', 'var')
            det.filename_pattern_scan = [p.detector.data_prefix];
            parts = strsplit(fill, del);
            if exist('burst', 'var')
                parts{burst} = '*';
            end
            if exist('pos', 'var')
                parts{pos} = '*';
            end
            det.filename_pattern_scan = [det.filename_pattern_scan strjoin(parts, del) det.file_extension];
        end
        
        input_vars = {};
        if isfield(det, 'filename_pattern_pos')
            k = 1;
            for jj=1:length(det.filename_pattern)
                switch det.filename_pattern{jj}.content
                    case 'pos'
                        continue
                    case 'burst'
                        input_vars{k} = det.filename_pattern{jj}.start;
                    case 'scan'
                        input_vars{k} = p.scan_number(p.scanID);
                end
                k = k+1;
            end
            filename_pattern_pos = fullfile(read_path, det.filename_pattern_pos);
            filename_pos = sprintf(filename_pattern_pos, input_vars{:});
            [~, pos_files] = find_files(filename_pos);
            numpos = size(pos_files,2);
            
            % apply natural sorting order i.e. sort 1,2,3,10,200 and not 1 10 100 2 20 200
            % important if the file makes are not defined as S%05i but rather S%i
            [~,idx] = natsort({pos_files.name});
            pos_files = pos_files(idx);
            
        end
        
        if isfield(det, 'filename_pattern_burst')
            k = 1;
            for jj=1:length(det.filename_pattern)
                switch det.filename_pattern{jj}.content
                    case 'pos'
                        input_vars{k} = det.filename_pattern{jj}.start;
                    case 'burst'
                        continue
                    case 'scan'
                        input_vars{k} = p.scan_number(p.scanID);
                end
                k = k+1;
                
            end
            filename_pattern_burst = fullfile(read_path, det.filename_pattern_burst);
            filename_burst = sprintf(filename_pattern_burst, input_vars{:});
            [~, burst_files] = find_files(filename_burst);
            numburst = size(burst_files,2);
        else
            numburst = 1;
        end
        
        
        
        if numburst > 1
            % if burst frames exist, we need to make sure that the file order is correct
            file_args = '[';
            for ii=1:length(det.filename_pattern)
                switch ii
                    case burst
                        file_args = [file_args ' det.filename_pattern{ii}.start + burstID-1'];
                    case pos
                        file_args = [file_args ' det.filename_pattern{ii}.start + posID-1'];
                    case scan
                        file_args = [file_args ' p.scan_number(p.scanID)'];
                end
            end
            file_args = [file_args ']'];
            for posID=1:numpos
                for burstID=1:numburst
                    files(burstID+(posID-1)*numburst).name = sprintf(det.filename_pattern_full, eval(file_args));
                end
            end
            datadir = read_path;
            
        else
            % if there are no burst frames, use the pos files
            datadir = read_path;
            files = pos_files;
        end
        
        if numel(files)==0
            error('Could not find any files using the filename pattern %s. \n ', fullfile(read_path, strrep(det.filename_pattern_full, '%', '%%')))
        end
    else
        % use wildcards
        files = find_files(fullfile(read_path, [det.filename_pattern det.file_extension]));
        if numel(files)==0
            error('Could not find any files using the filename pattern %s. \n ', fullfile(read_path, strrep(det.filename_pattern_full, '%', '%%')))
        end
    end
else
    % if no filename pattern was specified, just load everything containing the specified file extension
    [datadir, files] = find_files([read_path '*.' det.file_extension]);
    
    if numel(files)==0
        error('Could not find any files using the filename pattern %s.\n ', fullfile(read_path, ['*.' det.file_extension]))
    end
end


detStorage.files = [];

for ii=1:length(files)
    detStorage.files{ii} = fullfile(datadir, files(ii).name);
end

for ii=1:length(det.image_read_extraargs)
    if strcmpi(det.image_read_extraargs{ii}, 'H5Location')
        detStorage.h5_group{1} = det.image_read_extraargs{ii+1};
        break;
    end
end
%}
end

