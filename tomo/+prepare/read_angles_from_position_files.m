% Normally the tomography angles can be read from an angles file from OMNY,
% flOMNI, or LaMNI. However in one case this did not work and we had to
% extract the angle from the header of the ptycho positions file. This is a
% script that allows just that.
% Usage example:
% strpatt = '~/Data10/specES1/scan_positions/scan_%05d.dat';
% scannums = [271:277];
% angles = prepare.read_angles_from_position_files(strpatt,scannums)

function angles = read_angles_from_position_files(strpatt,scannums)

filenames = cell(numel(scannums),1);
angles = nan(numel(scannums),1);
for ii = 1:numel(scannums)
    filenames{ii} = sprintf(strpatt,scannums(ii));
    try
        aux = beamline.read_omny_pos(filenames{ii});
        angles(ii) = aux.lsamrot_encoder;
    catch
       fprintf('Scan %i is missing \n', scannums(ii)) 
    end
end
