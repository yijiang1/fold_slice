%CREATE_COLORMAP Create a colormap based on interpolation between anchor 
%values and colors
%
% cmap = plotting.create_colormap(values,colors,N)
%
% ** values    Array with numbers between 0 and 1
% ** colors    Array with [r g b] for each number in values
%
% *optional*
% ** N         Size of output colormap, default is 64
%
% returns
% ++ cmap      Nx3 array with the colormap
%
%
% see also: plotting.imagesc3D
%
% EXAMPLES:
%   A red colormap with 2 anchor colors, black and red
%   colormap(plotting.create_colormap([0 1],[0 0 0; 1 0 0]));

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2018 by Paul Scherrer Institute (http://www.psi.ch)    |
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

function cmap = create_colormap(values,colors,N)
if ~exist('N','var')
    N = 64;
end

v = linspace(0,1,N);
cmap = zeros(N,3);
[X V] = meshgrid([1 2 3],v);
cmap = interp2(repmat([1 2 3],[numel(values) 1]),repmat(values(:),[1 3]),colors,X,V);

end
