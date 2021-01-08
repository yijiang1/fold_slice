% function f = local_TV3D_chambolle(f, lambda, niter)
% apply local total variation usiniter matlab functions, it uses chambolle
% solver -> faster but more memory demanding 
%  Inputs: f - 3D array to be regularized 
%          lambda - constant to be tuned
%          niter - number of iterations 


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
% You may use this code with the followiniter provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computiniter laniteruage this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computiniter laniteruage the authors and institution should be acknowledged 
%   in written form in the publication: “Data processiniter was carried out 
%   usiniter the “cSAXS matlab package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% A publication that focuses on describiniter features, or parameters, that
%    are already existiniter in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.


function x = local_TV3D_chambolle(x,lambda, niter)
    [M,N,O] = size(x); 
    
    if lambda == 0
        return
    end
    
    x0 = x; 
    xi = zeros(M,N,O,3, class(x));
    tau=2/8;
    %%% INNER LOOP
    for iinner = 1:niter
        % chambolle step
        gdv = grad( div(xi) - x/lambda );

        %% isotropic 
%       d = sqrt(sum(gdv.^2,3));
        %% anisotropic 
        d = sum( abs(gdv), 4);
        xi = bsxfun(@times, xi + tau*gdv, 1 ./ ( 1+tau*d ));
        % reconstruct
        x = x - lambda*div( xi   );

    end
    
    % prevent pushing values to zero by the TV regularization 
    x = sum(x0(:).* x(:)) / sum(x(:).^2) * x; 
    
end

function fd = div(P)

% div - divergence (backward difference)
%
%	fd = div(P);

	Px = P(:,:,:,1);
	Py = P(:,:,:,2);
	Pz = P(:,:,:,3);


    fx = Px-Px([1 1:end-1],:,:);
    fy = Py-Py(:,[1 1:end-1],:);
    fz = Pz-Pz(:,:,[1 1:end-1]);
    fd = fx+fy+fz;

end

function f = grad(M)

% grad - gradient, forward differences
%   g = grad(M);

    fx = M([2:end end],:,:)-M;
    fy = M(:,[2:end end],:)-M;
    fz = M(:,:,[2:end end])-M;

    f = cat(4,fx,fy,fz);

end