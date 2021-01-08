function [U,S,V] = fsvd(A, k, i, usePowerMethod)
% FSVD Fast Singular Value Decomposition 
% 
%   [U,S,V] = FSVD(A,k,i,usePowerMethod) computes the truncated singular
%   value decomposition of the input matrix A upto rank k using i levels of
%   Krylov method as given in [1], p. 3.
% 
%   If usePowerMethod is given as true, then only exponent i is used (i.e.
%   as power method). See [2] p.9, Randomized PCA algorithm for details.
% 
%   [1] Halko, N., Martinsson, P. G., Shkolnisky, Y., & Tygert, M. (2010).
%   An algorithm for the principal component analysis of large data sets.
%   Arxiv preprint arXiv:1007.5510, 0526. Retrieved April 1, 2011, from
%   http://arxiv.org/abs/1007.5510. 
%   
%   [2] Halko, N., Martinsson, P. G., & Tropp, J. A. (2009). Finding
%   structure with randomness: Probabilistic algorithms for constructing
%   approximate matrix decompositions. Arxiv preprint arXiv:0909.4061.
%   Retrieved April 1, 2011, from http://arxiv.org/abs/0909.4061.
% 
%   See also SVD.
% 
%   Copyright 2011 Ismail Ari, http://ismailari.com.

    
    isSparse =  issparse(A);


    if nargin < 3
        i = 1;
    end

    % Take (conjugate) transpose if necessary. It makes H smaller thus
    % leading the computations to be faster
    if size(A,1) < size(A,2)
        A = A';
        isTransposed = true;
    else
        isTransposed = false;
    end

    n = size(A,2);
    extra_margin = 3;  % slighly improve precision 
    l = k + extra_margin;

    % Form a real n×l matrix G whose entries are iid Gaussian r.v.s of zero
    % mean and unit variance
    G = randn(n,l, 'single');


    if nargin >= 4 && usePowerMethod
        % Use only the given exponent
        H = A*G;
        for j = 2:i+1
            H = A * (A'*H);
        end
    else
        % Compute the m×l matrices H^{(0)}, ..., H^{(i)}
        % Note that this is done implicitly in each iteration below.
        
            if isSparse
                H = sparse(size(A,1), l * (i+1) );
            else
                H = zeros(size(A,1), l * (i+1), 'like', A);
            end
                        
            H(:,1:l) = A*G;

            for j = 2:i+1
               H(:,(j-1)*l +  (1:l)) = A * (A'*H(: , (j-2)*l +  (1:l)));
            end

        % Form the m×((i+1)l) matrix H
       
    end

    % Using the pivoted QR-decomposiion, form a real m×((i+1)l) matrix Q
    % whose columns are orthonormal, s.t. there exists a real
    % ((i+1)l)×((i+1)l) matrix R for which H = QR.  
    % XXX: Buradaki column pivoting ile yapılmayan hali.
    [Q,~] = qr(H,0);

   
    % Compute the n×((i+1)l) product matrix T = A^T Q
    T = A'*Q;

    % Form an SVD of T
    [Vt, St, W] = svd(T,'econ');

    % Compute the m×((i+1)l) product matrix
    Ut = Q*W;

    % Retrieve the leftmost m×k block U of Ut, the leftmost n×k block V of
    % Vt, and the leftmost uppermost k×k block S of St. The product U S V^T
    % then approxiamtes A. 

    
    if isTransposed
        V = Ut(:,1:k);
        U = Vt(:,1:k);     
    else
        U = Ut(:,1:k);
        V = Vt(:,1:k);
    end
    S = single(St(1:k,1:k));
end



