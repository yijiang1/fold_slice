% Generates a 1D orthonormal polynomial base
% polys = legendrepoly1D_2(X,maxorder,w);
% The weighting function has not been tested extensively
% Manuel Guizar - March 10, 2009

% Copyright (c) 2016, Manuel Guizar Sicairos, James R. Fienup, University of Rochester
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the University of Rochester nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function polys = legendrepoly1D_2(X,maxorder,w)

if nargin < 3
    w = 1;
end

[nc nr] = size(X);

polys = ones(nc,nr);
%%% Generation of polynomials
for ii = 0:maxorder,
    
    polys(:,:,ii+1) = (X.^(ii));   

end

%%% Normalization
for ii = 1:length(polys(1,1,:)),
%     polys(:,:,ii) = polys(:,:,ii)/sqrt(sum(sum(abs(polys(:,:,ii)).^2)));
    polys(:,:,ii) = polys(:,:,ii)/sqrt(sum(sum(w.*abs(polys(:,:,ii)).^2)));
end

%%% Orthonormalization
for ii = 2:length(polys(1,1,:)),
    for jj = 1:ii-1,
        polys(:,:,ii) = polys(:,:,ii) - sum(sum(polys(:,:,ii).*polys(:,:,jj).*w))*polys(:,:,jj);
    end
     
    polys(:,:,ii) = polys(:,:,ii)/sqrt(sum(sum(w.*polys(:,:,ii).^2)));
end
% polys(:,:,3) = polys(:,:,3) - sum(sum(polys(:,:,3).*polys(:,:,1)))*polys(:,:,1);
% polys(:,:,3) = polys(:,:,3)/sum(sum(polys(:,:,3).^2));

