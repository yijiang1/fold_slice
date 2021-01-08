%hdf5_pos_aps loads APS' data positions from hdf5 files (generated from python
%script)
%Written by YJ

function [ p ] = hdf5_pos_aps( p )

for ii = 1:p.numscans
    positions_real = zeros(0,2); 

    switch p.scan.type
        case 'custom'
            if isempty(p.scan.custom_positions_source) %guess the position file name from base path
                pos_file = strcat(p.base_path,sprintf(p.scan.format, p.scan_number(ii)),'/data_roi',p.scan.roi_label,'_para.hdf5');   
            else
                pos_file = p.scan.custom_positions_source;
            end
            if exist(pos_file,'file')
                
            	ppX = h5read(pos_file,'/ppX');
                ppY = h5read(pos_file,'/ppY');
                ppX = ppX(:);
                ppY = ppY(:);
                positions_real = zeros(length(ppX),2); 

                positions_real(:,1) = -ppY;
                positions_real(:,2) = -ppX;                
            else
                disp(strcat(p.base_path,sprintf(p.scan.format, p.scan_number(ii)),'/data_roi',p.scan.roi_label,'_para.hdf5'))
            	error('Could not find function or data file %s', pos_file);
            end
        case 'pre_recon'
            if ~isempty(p.scan.custom_positions_source) %guess the position file name from base path
                pos_file = p.scan.custom_positions_source;
            else
            	error('Position file is not given');
            end
            
            try
                r_output = load(pos_file,'outputs');
                r_p = load(pos_file,'p');
            	ppX = r_output.outputs.probe_positions(:,1)*r_p.p.dx_spec(1);
                ppY = r_output.outputs.probe_positions(:,2)*r_p.p.dx_spec(2);
                ppX = ppX(:);
                ppY = ppY(:);
                positions_real = zeros(length(ppX),2); 

                positions_real(:,1) = -ppY;
                positions_real(:,2) = -ppX;                
            catch
                error('Failed to load positions from %s', pos_file);
            end
        otherwise
            error('Unknown scan type %s.', p.scan.type);
    end
    utils.verbose(2, strcat('Loaded scan positions from:', pos_file))
    %scatter(ppX,ppY,'.');
    p.numpts(ii) = size(positions_real,1);
    p.positions_real = [p.positions_real ; positions_real]; %append position
end
    
end

