% PROP_FREE_NF    Near field propagation
%
%    wout = prop_free_nf(win, lambda, z, pixsize) returns the propagated wavefield
%    WIN by a distance Z, using wavelength LAMBDA. PIXSIZE is the dimension
%    of one pixel.
%
%    wout = prop_free_nf(win, lambda, z) is the same as above, assuming PIXSIZE=1
%    (that is, Z and LAMBDA are expressed in pixel units).
%
%   [~,H] = prop_free_nf(win, lambda, z, pixsize)    Does not compute propagation only 
%   computes and returns H, the Fourier domain propagation function
    
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

function [wout, H] = prop_free_nf(win, lambda, z, dx)

    if ndims(win) < 2
        error('Input wavefield should be at least 2-dimensional array!')
    end
    
    sz = size(win);

    if all(z==0)
        wout = win; 
        H = ones(sz(1),sz(2),'like',win);
        return
    end

    
    if nargin < 4
        dx = [1. 1.];
    end

    
    % Evaluate if aliasing could be a problem

    if any(sqrt(sum(1./dx.^2) * sum(1./(sz(1:2).*dx).^2)) * abs(z) * lambda > 1)
        utils.verbose(0,'Warning: there could be some aliasing issues...');
        utils.verbose(0,'(you could enlarge your array, or try a far field method)');
    end
    
    type = ones(1,'like',win);
    xgrid = type*[-sz(1)/2:floor((sz(1)-1)/2)]/sz(1)*lambda/dx(1); 
    ygrid = type*[-sz(2)/2:floor((sz(2)-1)/2)]/sz(2)*lambda/dx(min(end,2));
    [x1,x2] = ndgrid(xgrid, ygrid);
    q2 = math.fftshift_2D(x1.^2 + x2.^2);  % prevents numerical errors by keeping q2 unitless 
    
    z = reshape(z,1,1, []); % if more than one "z" is used, propagate each layer differenly
    
    if isa(win,'double')
        H = exp(2i * pi * (z / lambda) .* (sqrt(1 - q2) - 1)); 
    else
        % for single precision use approximation to avoid numerical issues,
        % precise with relative error less than 1e-7 for common use
        H = exp((2i * pi * (z / lambda)) .* (-q2/2)); 
    end
    
    if nargout < 2
        wout = ifft2(fft2(win) .* H  );
    else
        wout = [];  % do not calculate the wout and return only H factor 
    end
end

