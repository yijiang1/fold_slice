% Proyects a 1D function onto orthonormalized base, returns residual too
% The weighting function has not been tested extensively
% [coeffs reconstrproj] = projectleg1D_2(input,maxorder,Xext,w);
% March 10, 2009

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

function [coeffs, reconstrproj] = projectleg1D_2(input,maxorder,Xext,w)
import math.legendrepoly1D_2

polys = legendrepoly1D_2(Xext,maxorder,w);
reconstrproj = input; 


  
 
for ii = 1:length(polys(1,1,:)), 
    coeffs(ii) = sum(sum(reconstrproj.*polys(:,:,ii).*w));
% end
% 
% 
% for ii = 1:length(polys(1,1,:)), 
    reconstrproj = reconstrproj-polys(:,:,ii)*coeffs(ii);

end
 
% coeffs,

% figure(5);
% imagesc((real(reconstrproj)));
% axis square;
% colorbar;
% colormap gray;
% title('extracted from autocorrelation and projected to legendres')
% 
% figure(6);
% % imagesc((real(proj)));
% imagesc((real(reconstrproj-)));
% axis square;
% colorbar;
% colormap gray;
% title('projection to legendres')