function [pos_x, pos_y, mass, mu, sigma] = center(X, use_shift)
% -----------------------------------------------------------------------
% This file is part of the PTYCHOMAT Toolbox
% Author: Michal Odstrcil, 2016
% License: Open Source under GPLv3
% Contact: ptychomat@gmail.com
% Website: https://bitbucket.org/michalodstrcil/ptychomat
% -----------------------------------------------------------------------
% Description:  find center of mass of matrix X,  calculate variance if
% needed 
% inputs:
%   X  2D stacked images 
%   use_shift, if true, CoM will be calculated relatively to the center
%   of the image, default == true 



    if nargin < 2
        use_shift = true;
    end

    [N,M,~] = size(X);
    mass = sum(sum(X));
    xgrid = (1:M);
    ygrid = (1:N)';
    pos_x =  (sum(sum(X,1).*xgrid)./mass);
    pos_y =  (sum(sum(X,2).*ygrid)./mass);
    if nargout > 3
        mu = [pos_x, pos_y]';
    end
    
    if nargout == 5
        pos_xx = (sum(sum(X,1).*((1:M)-pos_x).^2)/mass);
        pos_yy = (sum(sum(X,2).*((1:N)-pos_y)'.^2)/mass);
        pos_xy = ((X*((1:M)-pos_x)')'*((1:N)-pos_y)'/mass);
        pos_yx = ((X*((1:M)-pos_x)')'*((1:N)-pos_y)'/mass);
        sigma =  [ pos_xx, pos_xy; pos_yx, pos_yy];
    end
    

    if use_shift
        pos_x = pos_x - M/2-0.5;  % defined so that center(ones(N),true) == [0,0]
        pos_y = pos_y - N/2-0.5;
    end
              
end
