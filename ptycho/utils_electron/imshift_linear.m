% IMSHIFT_LINEAR  will apply shift that can be different for 
%     each frame.  
%     + compared to imshift_fft, it does not have periodic boundary 
%     + it is based on linear interpolation, so it can be run fast on GPU 
%     + integer shift is equivalent to imshift_fft (up to the boundary condition)
%     - it needs for-loop for each frame -> it gets slow on GPU for
%       shifting my small images. In that case imshift_fft can be faster. 
%
%   img = imshift_linear(img, x,y, method)
%
%   Inputs:
%       **img   input image / stack of images 
%       **x     applied shift or vector of shifts for each frame 
%       **y     applied shift or vector of shifts for each frame 
%       **method choose interpolation method: nearest, {linear}, cubic , circ
% 
%   *returns*: 
%       ++img  shifted image / stack of images 
% 
%  see also: utils.imshift_fast, utils.imshift_fft




%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
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



    
function img = imshift_linear(img, x,y, method)


    if nargin  < 3
        y = x(:,2);
        x = x(:,1);
    end
    if nargin < 4
        method = 'linear';
    end
    
    if all(x==0) && all(y==0)
        return
    end

    
    real_img = isreal(img);
    [Nx, Ny,Nlayers] = size(img);

    if isscalar(x)
        x = ones(Nlayers,1) * x; 
    end
    if isscalar(y)
        y = ones(Nlayers,1) * y; 
    end
    
    if strcmpi(method, 'circ')
        % perform fast shift with circular boundary condition 
        X = 1:Nx;
        Y = 1:Ny;
        for ii = 1:Nlayers
            img(:,:,ii) = img(circshift(X,round(y(ii))), ...
                                              circshift(Y,round(x(ii))),ii);
        end
    else    

        for ii = 1:Nlayers
            %x(ii)
            %y(ii)
            img(:,:,ii) = interp2(single(img(:,:,ii)), single(-x(ii)+(1:Ny)),single(-y(ii)+(1:Nx)'), method,0);
        end
    
    end

    
  
  
end
