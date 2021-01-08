% mdaload.m  -  MATLAB/Octave routine for loading MDA files

% Copyright (c) 2016 UChicago Argonne, LLC,
%     as Operator of Argonne National Laboratory.
% This file is distributed subject to a Software License Agreement
% found in file LICENSE that is included with this distribution. 

%  Written by Dohn A. Arms, Argonne National Laboratory
%  Send comments to dohnarms@anl.gov

% 1.0 -- July 2016
%        Initial version

function mda = mdaload(filename)
[fileID,errmsg] = fopen( filename,'r','b');
if fileID < 0
    error( errmsg)
end
mda.version = float32grab(fileID);
ver = round(mda.version*100.0);
if ver ~= 140 && ver ~= 130 && ver ~= 120
    error('Incorrect MDA version')
end
mda.scan_number = int32grab(fileID);
mda.data_rank = int16grab(fileID);
if mda.data_rank < 1
    error('Internal inconsistency')
end
[mda.dimensions,cnt] = fread(fileID,mda.data_rank,'int32=>int32');
if (cnt ~= mda.data_rank) || any(mda.dimensions < 1)
    error('Internal inconsistency')
end
mda.regular = int16grab(fileID);
extra_offset = int32grab(fileID);
mda.scan = scangrab( fileID,mda.data_rank);
if extra_offset > 0
    fseek( fileID, extra_offset, 'bof');
    mda.extra = extragrab(fileID);
else
    mda.extra = [];
end
fclose(fileID);
end

function scan = scangrab(fileID,data_rank)
scan.scan_rank = int16grab(fileID);
if scan.scan_rank ~= data_rank
    error('Internal inconsistency')
end
scan.requested_points = int32grab(fileID);
if scan.requested_points < 1
    error('Internal inconsistency')
end
scan.last_point = int32grab(fileID);
if (scan.last_point < 0) || (scan.last_point > scan.requested_points)
    error('Internal inconsistency')
end
if scan.scan_rank > 1
    [offsets,cnt] = fread(fileID,scan.requested_points,'int32=>int32');
    if (cnt ~= scan.requested_points) || any(offsets(1:scan.last_point) == 0) || any(offsets < 0)
        error('Internal inconsistency')
    end
end
scan.name = strgrab(fileID);
scan.time = strgrab(fileID);
scan.number_positioners = int16grab(fileID);
if scan.number_positioners < 0
    error('Internal inconsistency')
end
scan.number_detectors = int16grab(fileID);
if scan.number_detectors < 0
    error('Internal inconsistency')
end
scan.number_triggers = int16grab(fileID);
if scan.number_triggers < 0
    error('Internal inconsistency')
end
if scan.number_positioners > 0
    for i = 1:scan.number_positioners
        scan.positioners(i) = posgrab(fileID);
    end
else
    scan.positioners = [];
end
if scan.number_detectors > 0
    for i = 1:scan.number_detectors
        scan.detectors(i) = detgrab(fileID);
    end
else
    scan.detectors = [];
end
if scan.number_triggers > 0
    for i = 1:scan.number_triggers
        scan.triggers(i) = triggrab(fileID);
    end
else
   scan.triggers = []; 
end
[scan.positioners_data,cnt] = fread(fileID, [scan.requested_points,scan.number_positioners],'float64=>float64');
if cnt ~= scan.requested_points*int32(scan.number_positioners)
    error('Internal inconsistency')
end
[scan.detectors_data,cnt] = fread(fileID, [scan.requested_points,scan.number_detectors],'float32=>float32');
if cnt ~= scan.requested_points*int32(scan.number_detectors)
    error('Internal inconsistency')
end
if scan.last_point  < scan.requested_points
    scan.positioners_data(scan.last_point+1:end,:) = [];
    scan.detectors_data(scan.last_point+1:end,:) = [];
end
if scan.scan_rank > 1
    for i = 1:scan.last_point
        if fseek( fileID, offsets(i), 'bof') == -1
            error('Premature end of file');
        end
        scan.sub_scans(i) = scangrab(fileID,scan.scan_rank-1);
    end
else
    scan.subscans = [];
end
end


function pos = posgrab(fileID)
pos.number = int16grab(fileID);
pos.name = strgrab(fileID);
pos.description = strgrab(fileID);
pos.step_mode = strgrab(fileID);
pos.unit = strgrab(fileID);
pos.readback_name = strgrab(fileID);
pos.readback_description = strgrab(fileID);
pos.readback_unit = strgrab(fileID);
end

function det = detgrab(fileID)
det.number = int16grab(fileID);
det.name = strgrab(fileID);
det.description = strgrab(fileID);
det.unit = strgrab(fileID);
end

function trig = triggrab(fileID)
trig.number = int16grab(fileID);
trig.name = strgrab(fileID);
trig.command = float32grab(fileID);
end

function extra = extragrab(fileID)
extra.number_pvs = int16grab(fileID);
    for i = 1:extra.number_pvs
        extra.pvs(i) = pvgrab(fileID);
    end
end

function pv = pvgrab(fileID)
pv.name = strgrab(fileID);
pv.descr = strgrab(fileID);
type = int16grab(fileID);
switch type
    case 0
        pv.type = 'char';
    case 29
        pv.type = 'int16';
    case 30
        pv.type = 'float32';
    case 32
        pv.type = 'int8';
    case 33
        pv.type = 'int32';
    case 34
        pv.type = 'float64';
    otherwise
        error('Internal inconsistency')
end
if type > 0
    pv.count = int16grab(fileID);
    if pv.count < 1
        error('Internal inconsistency')
    end
    pv.unit = strgrab(fileID);
else
    pv.count = 1;
    pv.unit = '';
end
switch type
    case 0
        pv.values = strgrab(fileID);
        pv.count = length(pv.values);
    case 29
        [pv.values,cnt] = fread(fileID,pv.count,'int16=>int16');
        if cnt ~= pv.count
            error('Internal inconsistency')
        end
        if mod(pv.count,2) == 1
            [~,cnt] = fread(fileID,1,'int16=>int16');
            if cnt ~= 1
                error('Internal inconsistency')
            end
        end
    case 30
        [pv.values,cnt] = fread(fileID,pv.count,'float32=>float32');
        if cnt ~= pv.count
            error('Internal inconsistency')
        end
    case 32
        [pv.values,cnt] = fread(fileID,pv.count,'int8=>int8');
        if cnt ~= pv.count
            error('Internal inconsistency')
        end
        len = mod(pv.count,4);
        if len > 0
            [~,cnt] = fread(fileID,4-len,'int8=>int8');
            if cnt ~= 4-len
                error('Internal inconsistency')
            end
        end
    case 33
        [pv.values,cnt] = fread(fileID,pv.count,'int32=>int32');
        if cnt ~= pv.count
            error('Internal inconsistency')
        end
    case 34
        [pv.values,cnt] = fread(fileID,pv.count,'float64=>float64');
        if cnt ~= pv.count
            error('Internal inconsistency')
        end
end
end

function val = int32grab(fileID)
[val,cnt] = fread(fileID,1,'int32=>int32');
if cnt ~= 1
    error('Premature end of file');
end
end

function val = int16grab(fileID)
[val,cnt] = fread(fileID,1,'int32=>int16');
if cnt ~= 1
    error('Premature end of file');
end
end

function val = float32grab(fileID)
[val,cnt] = fread(fileID,1,'float32=>float32');
if cnt ~= 1
    error('Premature end of file');
end
end

function str = strgrab(fileID)
[len1,cnt]=fread(fileID,1,'int32=>int16');
if cnt ~= 1
    error('Premature end of file');
end
if len1 > 0
    [len2,cnt]=fread(fileID,1,'int32=>int16');
    if cnt ~= 1
        error('Premature end of file');
    end
    [str,cnt]=fread(fileID,[1,len2],'char=>char');
    if cnt ~= len2
        error('Premature end of file');
    end
    m = mod(len2,4);
    if m > 0
        [~,cnt]=fread(fileID,4-m,'char=>char');
       if cnt ~= 4-m
           error('Premature end of file');
       end
    end
else
    str='';
end
end



