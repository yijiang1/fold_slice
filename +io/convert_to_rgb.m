function [rgb_data] = convert_to_rgb(data)
%Convert complex data into rgb image showing both magnitude and phase
%   Detailed explanation goes here

import math.sp_quantile

[W,H] = size(data);
adata = abs(data);

alpha  = 1e-3;
tmp= sort(adata(:));
MAX = tmp(ceil(end*(1-alpha)));
ind = adata > MAX;
data(ind) = MAX * data(ind) ./ abs(data(ind));
adata = abs(data);
range = sp_quantile(adata(:), [1e-2, 1-1e-2],10);
adata = (adata  - range(1) ) ./ ( range(2) - range(1) );

ang_data = angle(data); 
hue =  mod(ang_data+2.5*pi, 2*pi)/(2*pi);
hsv_data = [ hue(:) , ones(W*H,1), adata(:) ];

hsv_data = min(max(0, hsv_data),1);


rgb_data = hsv2rgb(hsv_data);

rgb_data = reshape(rgb_data, W,H,3);      
rgb_data = min(1,rgb_data);

end

