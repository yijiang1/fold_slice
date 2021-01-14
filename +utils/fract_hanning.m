% fract_hanning(outputdim,unmodsize)
% out = Square array containing a fractional separable Hanning window with
% DC in upper left corner.
% outputdim = size of the output array
% unmodsize = Size of the central array containing no modulation.
% Creates a square hanning window if unmodsize = 0 (or ommited), otherwise the output array 
% will contain an array of ones in the center and cosine modulation on the
% edges, the array of ones will have DC in upper left corner.

% February 8, 2007

% Slight update on August 17, 2009
% Added a warning

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

function out = fract_hanning(outputdim,unmodsize);

if nargin > 2,
    error('Too many input arguments'),
elseif nargin == 1,
    unmodsize = 0;
end

if outputdim < unmodsize,
    error('Output dimension must be smaller or equal to size of unmodulated window'),
end

if unmodsize<0,
    unmodsize = 0;
    warning('Specified unmodsize<0, setting unmodsize = 0')
end

N = [0:outputdim-1];
% N = ifftshift([-floor(outputdim/2):ceil(outputdim/2)-1]);
% N = [-floor(outputdim/2):ceil(outputdim/2)-1];
[Nc,Nr] = meshgrid(N,N);

if unmodsize == 0,
    out = (1+cos(2*pi*Nc/outputdim)).*(1+cos(2*pi*Nr/outputdim))/4;
else
    % Columns modulation
    out = (1+cos(2*pi*(Nc- floor((unmodsize-1)/2) )/(outputdim+1-unmodsize)))/2;
    if floor((unmodsize-1)/2)>0,
        out(:,1:floor((unmodsize-1)/2)) = 1;
    end
    out(:,floor((unmodsize-1)/2) + outputdim+3-unmodsize:length(N)) = 1;
    % Row modulation
    out2 = (1+cos(2*pi*(Nr- floor((unmodsize-1)/2) )/(outputdim+1-unmodsize)))/2;
    if floor((unmodsize-1)/2)>0,
        out2(1:floor((unmodsize-1)/2),:) = 1;
    end
    out2(floor((unmodsize-1)/2) + outputdim+3-unmodsize:length(N),:) = 1;
    
    out = out.*out2;
end  
% out = ifftshift(out);
% one-edge at Nc = floor((unmodsize-1)/2)
% other-edge at Nc = floor((unmodsize-1)/2) + (outputdim+1-unmodsize)
%%% FINISH UP THIS CODE TO DO THE LOW PASS RECONSTRUCTION

return;