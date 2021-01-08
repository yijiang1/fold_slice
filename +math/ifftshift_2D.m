function [x,idx]  = ifftshift_2D(x)
% -----------------------------------------------------------------------
% This file is part of the PTYCHOMAT Toolbox
% Author: Michal Odstrcil, 2016
% License: Open Source under GPLv3
% Contact: ptychomat@gmail.com
% Website: https://bitbucket.org/michalodstrcil/ptychomat
% -----------------------------------------------------------------------
% Description:   faster version of matlab fftshift to work for stack of
% 2D images 
% Inputs:  2D or stack of 2D images 
% Outputs: 
%   x - 2D or stack of 2D images after fftshift along first 2 dimensions
%   idx - precalculated indices for fftshift operation


numDims = 2;
idx = cell(1, numDims);
for k = 1:numDims
    m = size(x, k);
    p = floor(m/2);
    idx{k} = [p+1:m 1:p];
end
x = x(idx{:},:,:);


end
