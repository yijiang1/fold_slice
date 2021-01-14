% fract_hanning_pad(outputdim,filterdim,unmodsize)
% out = Square array containing a fractional separable Hanning window with
% DC in upper left corner.
% outputdim = size of the output array
% filterdim = size of filter (it will zero pad if filterdim<outputdim
% unmodsize = Size of the central array containing no modulation.
% Creates a square hanning window if unmodsize = 0 (or ommited), otherwise the output array 
% will contain an array of ones in the center and cosine modulation on the
% edges, the array of ones will have DC in upper left corner.

% August 17, 2009

% Copyright (c) 2016, Manuel Guizar Sicairos, University of Rochester
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

function out = fract_hanning_pad(outputdim,filterdim,unmodsize);
import utils.fract_hanning

if nargin > 3,
    error('Too many input arguments'),
elseif nargin == 1,
    unmodsize = 0;
    filterdim = outputdim;
end

if outputdim < unmodsize,
    error('Output dimension must be smaller or equal to size of unmodulated window'),
end

if outputdim < filterdim,
    error('Filter cannot be larger than output size'),
end

if unmodsize<0,
    unmodsize = 0;
    warning('Specified unmodsize<0, setting unmodsize = 0')
end

out = zeros(outputdim);
out(round(outputdim/2+1-filterdim/2):round(outputdim/2+1+filterdim/2-1),...
    round(outputdim/2+1-filterdim/2):round(outputdim/2+1+filterdim/2-1)) ...
    = fftshift(fract_hanning(filterdim,unmodsize));
out = fftshift(out);

return;