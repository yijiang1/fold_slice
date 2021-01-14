function [mask] = make_circular_mask(N, radius)
%Make a circular mask. Made by YJ
%   N: size of image
%   radius: radius
if length(N)==1
    x = linspace(-floor(N/2),ceil(N/2)-1,N);
    y = linspace(-floor(N/2),ceil(N/2)-1,N);
else
    x = linspace(-floor(N(2)/2),ceil(N(2)/2)-1,N(2));
    y = linspace(-floor(N(1)/2),ceil(N(1)/2)-1,N(1));
end
[Y, X] = meshgrid(x,y);
S = sqrt(X.^2+Y.^2);

mask = S <= radius;
end

