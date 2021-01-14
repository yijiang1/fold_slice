% filt2d_pad(outputdim,filterdim,unmodsize)
% out = Square array containing a fractional separable Hanning window with
% DC in upper left corner.
% outputdim = size of the output array
% filterdim = size of filter (it will zero pad if filterdim<outputdim
% unmodsize = Size of the central array containing no modulation.
% Creates a square hanning window if unmodsize = 0 (or ommited), otherwise the output array 
% will contain an array of ones in the center and cosine modulation on the
% edges, the array of ones will have DC in upper left corner.
% Code based in fract_hanning_pad

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
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing language this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS matlab package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.


% License for fract_hanning_pad:
% Manuel Guizar - August 17, 2009
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

function out = filt2d_pad(outputdim,filterdim,unmodsize, varargin)
import utils.filt2d

if nargin == 1
    unmodsize = 0;
    filterdim = outputdim;
end

if any(outputdim < unmodsize)
    error('Output dimension must be smaller or equal to size of unmodulated window'),
end

if any(outputdim < filterdim)
    error('Filter cannot be larger than output size'),
end

if any(unmodsize<0)
    unmodsize = [0 0];
    warning('Specified unmodsize<0, setting unmodsize = 0')
end

if length(outputdim)<2
    outputdim = [outputdim outputdim];
elseif length(outputdim)>2
    error('3D filters are not supported.')
end

if length(unmodsize)<2
    unmodsize = [unmodsize unmodsize];
elseif length(unmodsize)>2
    error('3D filters are not supported.')
end

if length(filterdim)<2
    filterdim = [filterdim filterdim];
elseif length(filterdim)>2
    error('3D filters are not supported.')
end

out = zeros(outputdim);
out(round(outputdim(1)/2+1-filterdim(1)/2):round(outputdim(1)/2+1+filterdim(1)/2-1),...
    round(outputdim(2)/2+1-filterdim(2)/2):round(outputdim(2)/2+1+filterdim(2)/2-1)) ...
    = fftshift(filt2d(filterdim,unmodsize,varargin{:}));
out = fftshift(out);

return;