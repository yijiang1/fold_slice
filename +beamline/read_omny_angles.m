% read_omny_angles( OMNY_angles_file, scannums, tomo_id )
% OMNY_angles_file  -  File with Scan number, angle target, angle readout
% scannums          -  Array of scan numbers
% tomo_id           -  integer or list of integers, only if the scannums is empty
%
% out               - Contains fields with scan, target_angle, readout_angle
% errorflag         - = 1 if at least one scan was not found

function [ out errorflag ] = read_omny_angles( OMNY_angles_file, scannums, tomo_id, scannum_offset )

if ~exist(OMNY_angles_file, 'file')
    error('Missing OMNY file: %s', OMNY_angles_file)
end

if ~exist('tomo_id')
    tomo_id = [];
end

if (~isempty(scannums))&&(~isempty(tomo_id))
    error('You have provided both scannums and tomo_id, please provide just either scannums OR tomo_id. One of them should be empty ( =[] ).')
end
if (isempty(scannums))&&(isempty(tomo_id))
    error('You have not provided scannums or tomo_id, please provide either scannums OR tomo_id. One of them should be empty ( =[] ).')
end
fid = fopen(OMNY_angles_file);

% check omny file type
ln = fgetl(fid);
switch numel(strsplit(ln, ' ')) 
    case {3,6}
        outmat = textscan(fid,'%f %f %f %f %f %s');
        fclose(fid);
        out = [];
        errorflag = 0;
        counter = 1;
        
        for ii = 1:numel(scannums)
            ind = find(outmat{1}==scannums(ii)+scannum_offset,1,'last');
            if isempty(ind)
                fprintf('Did not find Scan %d in %s\n',scannums(ii),OMNY_angles_file);
                errorflag = 1;
            else
                out.scan(counter) = outmat{1}(ind)-scannum_offset;
                out.target_angle(counter) = outmat{2}(ind);
                out.readout_angle(counter) = outmat{3}(ind);
                out.subtomo_num(counter) = outmat{4}(ind);
                out.detpos_num(counter) = outmat{5}(ind);
                out.sample_name(counter) = outmat{6}(ind);
                counter = counter+1;
            end
        end
    case 7
        outmat = textscan(fid,'%f %f %f %f %f %f %s');
        fclose(fid);
        out = [];
        errorflag = 0;
        counter = 1;
        if ~isempty(scannums)
            for ii = 1:numel(scannums)
                ind = find(outmat{1}==scannums(ii)+scannum_offset,1,'last');
                if isempty(ind)
                    fprintf('Did not find Scan %d in %s\n',scannums(ii),OMNY_angles_file);
                    errorflag = 1;
                else
                    out.scan(counter) = outmat{1}(ind)-scannum_offset;
                    out.target_angle(counter) = outmat{2}(ind);
                    out.readout_angle(counter) = outmat{3}(ind);
                    out.tomo_id(counter) = outmat{4}(ind);
                    out.subtomo_num(counter) = outmat{5}(ind);
                    out.detpos_num(counter) = outmat{6}(ind);
                    out.sample_name(counter) = outmat{7}(ind);
                    counter = counter+1;
                end
            end
        elseif ~isempty(tomo_id) %Note: Not used at APS
            ind = find(ismember(outmat{4},tomo_id));
            if isempty(ind)
                fprintf(['Did not find tomo_id ',repmat('%i ',1,length(tomo_id)),' in %s\n'],tomo_id,OMNY_angles_file);
                errorflag = 1;
            end
            out.scan            = outmat{1}(ind);
            out.target_angle    = outmat{2}(ind);
            out.readout_angle   = outmat{3}(ind);
            out.tomo_id         = outmat{4}(ind);
            out.subtomo_num     = outmat{5}(ind);
            out.detpos_num      = outmat{6}(ind);
            out.sample_name     = outmat{7}(ind);
        end
    otherwise
        error('Unknown OMNY file format.')
end


return
end

