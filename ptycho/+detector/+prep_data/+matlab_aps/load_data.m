%LOAD_DATA prepare filenames and load data
% ** p      p structure
%
% returns:
% ++ p      p structure
%
% see also: detector.prep_data.matlab_ps.prepare_data
% written by YJ


function [ p ] = load_data( p )
import utils.find_files
import io.image_read

det.params = p.detectors(p.scanID).params;
detStorage = p.detectors(p.scanID).detStorage;
%% prepare filenames
if isempty(detStorage.files)
    [p] = det.params.get_filename(p);
    if numel(detStorage.files) == 0
        error('Did not find any files.')
    end
end

utils.verbose(2, 'Loading raw data of scan %05d.', p.scan_number(p.scanID))

files = detStorage.files;
utils.verbose(2, strcat('HDF5 file name:', files{p.scanID}))
%disp(files{p.scanID})

try
    data = h5read(files{p.scanID},'/dp');
catch
    error('Failed to load dp from %s', files{p.scanID});
end
data = squeeze(data);
utils.verbose(2, strcat('Loaded data from:', files{p.scanID}))

if det.params.orientation(1)
    utils.verbose(2, 'Transposing diffraction patterns')
    data = permute(data, [2 1 3]);
end
if det.params.orientation(2) && det.params.orientation(3)
    utils.verbose(2, 'Flipping diffraction patterns')

    data = rot90(data,2);  % merge the fliplr and flipup operations 
else
    if det.params.orientation(2)
        data = fliplr(data);
    end
    if det.params.orientation(3)
        data = flipud(data);
    end
end

%{
%% select data loading routine
current_version = version;
ver_str = strsplit(current_version, '.');
ver_num = str2double([ver_str{1} '.' ver_str{2} ver_str{3}]);
if ver_num >= 9.4 
    c_reader = true;
else
    utils.verbose(3, 'Fast data reader is not available for Matlab version %d. Switching to image_read.', ver_num)
    c_reader = false;
end
if c_reader && ~ismember(det.params.file_extension, {'h5', 'cbf', 'tiff', 'tif'})
    utils.verbose(3, 'Fast data reader is not available for selected file format %s. Switching to image_read.', det.params.file_extension)
    c_reader = false;
end


%% load data
if c_reader
    % use fast data reader
    if strcmpi(det.params.file_extension, 'tif')
        det.params.file_extension = 'tiff'; 
    end
    if strcmpi(det.params.file_extension, 'tiff')
        % keep results consistent with matlab's imread and io.image_read for tiff files 
        det.params.orientation(1) = ~det.params.orientation(1); 
    end
    % convert PtychoShelves center to raw data center
    arg.ctr = detStorage.ctr-1;
    sz = det.params.geometry.sz;

    if det.params.orientation(3)
        arg.ctr(1) = round(sz(1) - arg.ctr(1));
    end
    if det.params.orientation(2)
        arg.ctr(2) = round(sz(2) - arg.ctr(2));
    end  
    if det.params.orientation(1)
        arg.ctr = fliplr(arg.ctr);
    end
    if iscolumn(arg.ctr)
        arg.ctr = arg.ctr';
    end
    
    % flip XY
    arg.ctr = fliplr(arg.ctr);
    
    % create structure for c_reader
    arg.data_path = files;
    arg.nthreads = p.io.data_nthreads;
    arg.precision = p.io.data_precision;
    arg.extension = det.params.file_extension;
    arg.asize = detStorage.read_size;
    arg.data_location = detStorage.h5_group;
    assert(all(arg.ctr>0),'Raw data center position has to be positive.')
    % load data and permute
    utils.verbose(2, 'Loading raw data of scan S%05d.', p.scan_number(p.scanID))
    data = io.read_measurement(arg);
    data = squeeze(data);
    size(data)

    if det.params.orientation(1)
        data = permute(data, [2 1 3]);
    end
    if det.params.orientation(2) && det.params.orientation(3)
        data = rot90(data,2);  % merge the fliplr and flipup operations 
    else
        if det.params.orientation(2)
            data = fliplr(data);
        end
        if det.params.orientation(3)
            data = flipud(data);
        end
    end
    
    
    
else
    if isempty(detStorage.h5_group)
        if numel(files)==1
            % one directory contains single file
            dataaux = image_read(files{1}, det.params.image_read_extraargs);
            data = dataaux.data;
        else
            % multiple files per directory
            dataaux = image_read(files, det.params.image_read_extraargs);
            data = dataaux.data;
        end
    else
        data = zeros([p.asize numel(detStorage.h5_group)]);
        if numel(files)==1
            if numel(detStorage.h5_group)==1
                numel(detStorage.h5_group)
                % one hdf5 file; one group
                dataaux = image_read(files{1}, det.params.image_read_extraargs{:}, 'H5Location', detStorage.h5_group{1});
                data = dataaux.data;
            else
                % one hdf5 file; multiple groups
                for ii=1:numel(detStorage.h5_group)
                    dataaux = image_read(files{1}, det.params.image_read_extraargs{:}, 'H5Location', detStorage.h5_group{ii});
                    data(:,:,ii) = dataaux.data;
                end
            end
        else
            if numel(detStorage.h5_group)==1
                % multiple files, single H5 group
                dataaux = image_read(files, det.params.image_read_extraargs{:}, 'H5Location', detStorage.h5_group{1});
                data = dataaux.data;
            else
                % multiple files, but different H5 group
                for ii=1:numel(detStorage.h5_group)
                    dataaux = image_read(files{ii}, det.params.image_read_extraargs{:}, 'H5Location', detStorage.h5_group{ii});
                    data(:,:,ii) = dataaux.data;
                end
            end
        end
    end
end
%}

detStorage.data = data;
end
