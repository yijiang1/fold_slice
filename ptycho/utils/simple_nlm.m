function [output]=simple_nlm(input,t,f,h1,h2,selfsim)
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 %  input   : image to be filtered
 %  t       : radius of search window
 %  f       : radius of similarity window
 %  h1,h2   : w(i,j) = exp(-||GaussFilter(h1) .* (p(i) - p(j))||_2^2/h2^2)
 %  selfsim : w(i,i) = selfsim, for all i
 %
 %  Note:
 %    if selfsim = 0, then w(i,i) = max_{j neq i} w(i,j), for all i
 %
 %  Author: Christian Desrosiers
 %  Date: 07-07-2015
 %
 %  Reimplementation of the Non-Local Means Filter by Jose Vicente Manjon-Herrera
 %
 %  For details see:
 %     A. Buades, B. Coll and J.M. Morel, "A non-local algorithm for image denoising"
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 [m, n]=size(input);
 pixels = input(:);
 
 s = m*n;
 
 psize = 2*f+1;
 nsize = 2*t+1;

 % Compute patches
 padInput = padarray(input,[f f],'symmetric'); 
 filter = fspecial('gaussian',psize,h1);
 patches = repmat(sqrt(filter(:))',[s 1]) .* im2col(padInput, [psize psize], 'sliding')';
 
 % Compute list of edges (pixel pairs within the same search window)
 indexes = reshape(1:s, m, n);
 padIndexes = padarray(indexes, [t t]);
 neighbors = im2col(padIndexes, [nsize, nsize], 'sliding');
 TT = repmat(1:s, [nsize^2 1]);
 edges = [TT(:) neighbors(:)];
 RR = find(TT(:) >= neighbors(:));
 edges(RR, :) = [];
 
 % Compute weight matrix (using weighted Euclidean distance)
 diff = patches(edges(:,1), :) - patches(edges(:,2), :);
 V = exp(-sum(diff.*diff,2)/h2^2); 
 W = sparse(edges(:,1), edges(:,2), V, s, s);
 
 % Make matrix symetric and set diagonal elements
 if selfsim > 0
    W = W + W' + selfsim*speye(s);
 else
     maxv = max(W,[],2);
     W = W + W' + spdiags(maxv, 0, s, s);
 end     
 
 % Normalize weights
 W = spdiags(1./sum(W,2), 0, s, s)*W;
 
 % Compute denoised image
 output = W*pixels;
 output = reshape(output, m , n);