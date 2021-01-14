% read_omny_angles( OMNY_angles_file, scannums, tomo_id )
% OMNY_angles_file  -  File with Scan number, angle target, angle readout
% scannums          -  Array of scan numbers
% tomo_id           -  integer or list of integers, only if the scannums is empty
%
% out               - Contains fields with scan, target_angle, readout_angle
% errorflag         - = 1 if at least one scan was not found

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

function [ out errorflag ] = read_omny_angles( OMNY_angles_file, scannums, tomo_id )

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
            ind = find(outmat{1}==scannums(ii),1,'last');
            if isempty(ind)
                fprintf('Did not find Scan %d in %s\n',scannums(ii),OMNY_angles_file);
                errorflag = 1;
            else
                out.scan(counter) = outmat{1}(ind);
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
                ind = find(outmat{1}==scannums(ii),1,'last');
                if isempty(ind)
                    fprintf('Did not find Scan %d in %s\n',scannums(ii),OMNY_angles_file);
                    errorflag = 1;
                else
                    out.scan(counter) = outmat{1}(ind);
                    out.target_angle(counter) = outmat{2}(ind);
                    out.readout_angle(counter) = outmat{3}(ind);
                    out.tomo_id(counter) = outmat{4}(ind);
                    out.subtomo_num(counter) = outmat{5}(ind);
                    out.detpos_num(counter) = outmat{6}(ind);
                    out.sample_name(counter) = outmat{7}(ind);
                    counter = counter+1;
                end
            end
        elseif ~isempty(tomo_id)
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

