%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unwrapping phase based on Ghiglia and Romero (1994) based on weighted and unweighted least-square method
% URL: https://doi.org/10.1364/JOSAA.11.000107
% Inputs:
%   * psi: wrapped phase from -pi to pi
%   * weight: weight of the phase (optional, default: all ones)
% Output:
%   * phi: unwrapped phase from the weighted (or unweighted) least-square phase unwrapping
% Author: Muhammad F. Kasim (University of Oxford, 2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function phi = phase_unwrap(psi, weight)
    if (nargin < 2) % unweighted phase unwrap
        % get the wrapped differences of the wrapped values
        dx = [zeros([size(psi,1),1]), wrapToPi(diff(psi, 1, 2)), zeros([size(psi,1),1])];
        dy = [zeros([1,size(psi,2)]); wrapToPi(diff(psi, 1, 1)); zeros([1,size(psi,2)])];
        rho = diff(dx, 1, 2) + diff(dy, 1, 1);
        
        % get the result by solving the poisson equation
        phi = solvePoisson(rho);
        
    else % weighted phase unwrap
        % check if the weight has the same size as psi
        if (~all(size(weight) == size(psi)))
            error('Argument error: Size of the weight must be the same as size of the wrapped phase');
        end
        
        % vector b in the paper (eq 15) is dx and dy
        dx = [wrapToPi(diff(psi, 1, 2)), zeros([size(psi,1),1])];
        dy = [wrapToPi(diff(psi, 1, 1)); zeros([1,size(psi,2)])];
        
        % multiply the vector b by weight square (W^T * W)
        WW = weight .* weight;
        WWdx = WW .* dx;
        WWdy = WW .* dy;
        
        % applying A^T to WWdx and WWdy is like obtaining rho in the unweighted case
        WWdx2 = [zeros([size(psi,1),1]), WWdx];
        WWdy2 = [zeros([1,size(psi,2)]); WWdy];
        rk = diff(WWdx2, 1, 2) + diff(WWdy2, 1, 1);
        normR0 = norm(rk(:));
        
        % start the iteration
        eps = 1e-6;
        k = 0;
        phi = zeros(size(psi));
        while (~all(rk == 0))
            zk = solvePoisson(rk);
            k = k + 1;
            if (k == 1) pk = zk;
            else 
                betak = sum(sum(rk .* zk)) / sum(sum(rkprev .* zkprev));
                pk = zk + betak * pk;
            end
            
            % save the current value as the previous values
            rkprev = rk;
            zkprev = zk;
            
            % perform one scalar and two vectors update
            Qpk = applyQ(pk, WW);
            alphak = sum(sum(rk .* zk)) / sum(sum(pk .* Qpk));
            phi = phi + alphak * pk;
            rk = rk - alphak * Qpk;
            
            % check the stopping conditions
            if ((k >= numel(psi)) || (norm(rk(:)) < eps * normR0)) break; end;
        end
    end
end

function phi = solvePoisson(rho)
    % solve the poisson equation using dct
    dctRho = dct2(rho);
    [N, M] = size(rho);
    [I, J] = meshgrid([0:M-1], [0:N-1]);
    dctPhi = dctRho ./ 2 ./ (cos(pi*I/M) + cos(pi*J/N) - 2);
    dctPhi(1,1) = 0; % handling the inf/nan value
    
    % now invert to get the result
    phi = idct2(dctPhi);
    
end

% apply the transformation (A^T)(W^T)(W)(A) to 2D matrix
function Qp = applyQ(p, WW)
    % apply (A)
    dx = [diff(p, 1, 2), zeros([size(p,1),1])];
    dy = [diff(p, 1, 1); zeros([1,size(p,2)])];
    
    % apply (W^T)(W)
    WWdx = WW .* dx;
    WWdy = WW .* dy;
    
    % apply (A^T)
    WWdx2 = [zeros([size(p,1),1]), WWdx];
    WWdy2 = [zeros([1,size(p,2)]); WWdy];
    Qp = diff(WWdx2,1,2) + diff(WWdy2,1,1);
end
