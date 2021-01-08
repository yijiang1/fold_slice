%convert_APS_ptycho_recon.m

flip_probe = false;
%% load probe
probe_temp = outputs.probe{1};
probe = zeros(size(probe_temp,1),size(probe_temp,2),length(outputs.probe));
for i=1:length(outputs.probe)
   probe(:,:,i) = outputs.probe{i}(:,:,1,1);
   if flip_probe
       probe(:,:,i) = rot90(probe(:,:,i),2);    
   end
end
clear i probe_temp
%% add parameters
p = {};
p.binning = false;
p.detector.binning = false;