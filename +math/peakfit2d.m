% Find sub-sample location of a global peak within 2D-matrix by applying 
% two dimensional polynomial fit & extremum detection. 
%
% Sample usage: 
% >> M = exp(-((1:30) - 19.5).^2/(2*5^2)); % gauss: center=19.5; sigma=5
% >> P = peakfit2d(M'*M);                  % find peak in 2D-gauss
% >> disp(P);
%   19.5050   19.5050
%
% Algebraic solution derived with the following steps:
%
% 0.) Define Approximation-Function: 
%
%     F(x,y) => z = a*x^2+b*x*y+c*x+d+e*y^2+f*y
%
% 1.) Formulate equation for sum of squared differences with
%
%     x=-1:1,y=-1:1,z=Z(x,y)
%
%     SSD = [ a*(-1)^2+b*(-1)*(-1)+c*(-1)+d+e*(-1)^2+f*(-1) - Z(-1,-1) ]^2 + ...
%              ...
%             a*(+1)^2+b*(+1)*(+1)+c*(+1)+d+e*(+1)^2+f*(+1) - Z(-1,-1) ]^2
%        
% 2.) Differentiate SSD towards each parameter
%
%     dSSD / da = ...
%              ...
%     dSSD / df = ...
%
% 3.) Solve linear system to get [a..f]
%
% 4.) Differentiate F towards x and y and solve linear system for x & y
%
%     dF(x,y) / dx = a*... = 0 !
%     dF(x,y) / dy = b*... = 0 !

% Copyright (c) 2010, Eric
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
% * Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in
% the documentation and/or other materials provided with the distribution
% * Neither the name of the HTWK Leipzig nor the names
% of its contributors may be used to endorse or promote products derived
% from this software without specific prior written permission.
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
% POSSIBILITY OF SUCH DAMAGE

function P = peakfit2d(Z)
import math.peakfit2d

%% Check input
sZ = size(Z);
if min(sZ)<2
    disp('Wrong matrix size. Input matrix should be numerical MxN type.');
    P = [0 0];
    return;
end


%% peak approximation using 2D polynomial fit within 9 point neighbourship

% find global maximum and extract 9-point neighbourship
[v,p] = max(Z(:));
[yp,xp]=ind2sub(sZ,p); 
if (yp==1)||(yp==sZ(1))||(xp==1)||(xp==sZ(2))
    disp('Maximum position at matrix border. No subsample approximation possible.');
    P = [yp xp];
    return;
end
K = Z(yp-1:yp+1,xp-1:xp+1);

% approximate polynomial parameter
a = (K(2,1)+K(1,1)-2*K(1,2)+K(1,3)-2*K(3,2)-2*K(2,2)+K(2,3)+K(3,1)+K(3,3));
b = (K(3,3)+K(1,1)-K(1,3)-K(3,1));
c = (-K(1,1)+K(1,3)-K(2,1)+K(2,3)-K(3,1)+K(3,3));
%d = (2*K(2,1)-K(1,1)+2*K(1,2)-K(1,3)+2*K(3,2)+5*K(2,2)+2*K(2,3)-K(3,1)-K(3,3));
e = (-2*K(2,1)+K(1,1)+K(1,2)+K(1,3)+K(3,2)-2*K(2,2)-2*K(2,3)+K(3,1)+K(3,3));
f = (-K(1,1)-K(1,2)-K(1,3)+K(3,1)+K(3,2)+K(3,3));

% (ys,xs) is subpixel shift of peak location relative to point (2,2)
ys = (6*b*c-8*a*f)/(16*e*a-9*b^2);
xs = (6*b*f-8*e*c)/(16*e*a-9*b^2);

P = [ys+yp xs+xp];