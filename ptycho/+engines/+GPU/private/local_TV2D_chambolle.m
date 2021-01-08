% function f = local_TV2D_chambolle(f, lambda, niter)
% apply local total variation usiniter matlab functions, it uses chambolle
% solver -> faster but more memory demanding 
%  Inputs: x - 2D array to be regularized 
%          lambda - constant to be tuned
%          niter - number of iterations 
% Modified from PSI's tomo code. Written by Jonathan Schwartz at U. Mich

function x = local_TV2D_chambolle(x,lambda, niter)
    [M,N] = size(x); 
    
    if lambda == 0
        return
    end
    
    x0 = x; 
    
    xi = zeros(M,N,2, class(x));
    tau=1/8;

    %%% INNER LOOP
    for iinner = 1:niter
        % chambolle step
        gdv = grad( div(xi) - x/lambda );

        %% isotropic 
        d = sqrt(sum(gdv.^2,3));
        
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
%   fd = div(P);

    Px = P(:,:,1);
    Py = P(:,:,2);

    fx = Px-Px([1 1:end-1],:);
    fy = Py-Py(:,[1 1:end-1]);

    fd = fx+fy;

end

function f = grad(M)

% grad - gradient, forward differences
%   g = grad(M);

    fx = M([2:end end],:)-M;
    fy = M(:,[2:end end])-M;

    f = cat(3,fx,fy);

end