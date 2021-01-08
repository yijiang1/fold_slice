function ellipse_fit_positions(samposx,samposy,scans,lamni_fit_file)

% addpath ../
% 
% cf{1} = load('Claire_click_squarymcsquareface.mat');
% samposx_all = [];
% samposy_all = [];
% scans_all = [];
% for ii = 1:numel(cf)
%     samposx_all = [samposx_all cf{ii}.samposx];
%     samposy_all = [samposy_all cf{ii}.samposy];
%     scans_all = [scans_all cf{ii}.scans];
% end

figure(1);
clf
plot(samposx,samposy,'o')

%%
if exist('lamni_fit_file')
    [~,lamni_name,~] = fileparts(lamni_fit_file);
    %theta = prepare.read_angles_from_position_files('~/specES1/scan_positions/scan_%05d.dat',scans);
    theta = prepare.read_angles_from_position_files('/mnt/micdata2/lamni/2020-2/comm_33IDD_2/specES1/scan_positions/scan_%05d.dat',scans);
    
    theta_interp = [min(theta)-1:max(theta)+1];
    
    samposx_interp = interp1(theta,samposx,theta_interp,'spline');
    samposy_interp = interp1(theta,samposy,theta_interp,'spline');
    hold on
    plot(samposx_interp, samposy_interp,'--')
    hold off
    
    corr_filename = sprintf('correction_lamni_um_S%05d_%s.txt',scans(1),lamni_name);
    h = fopen(corr_filename,'w');
    fprintf(h,'corr_elements = %d \n',length(theta_interp));
    for ii = 1:length(samposx_interp)
        fprintf(h,'%s[%d] = %.6f \n','corr_angle',ii-1,theta_interp(ii));
        fprintf(h,'%s[%d] = %.6f \n','corr_pos_x',ii-1, samposx_interp(ii));
        fprintf(h,'%s[%d] = %.6f \n','corr_pos_y',ii-1, samposy_interp(ii));
    end
    fclose(h);
    fprintf('Wrote succesfully to %s\n',corr_filename)
end
