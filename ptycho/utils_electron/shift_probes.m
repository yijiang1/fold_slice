%shift_probes.m
sy = -21;
sx = 7;
probe_temp = zeros(128,128,size(probe,3),1);
for i=1:size(probe,3)
    probe_temp(:,:,i,1) = imshift_linear(probe(:,:,i,1), sx,sy);
    
    
end
    