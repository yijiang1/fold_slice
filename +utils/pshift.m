% AOUT = pshift(AIN, CTRPOS)
%
% Shift array AIN periodically so that CTRPOS is placed at (1,1). 
%

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

function aout = pshift(ain, ctrpos, varargin)

sz = size(ain);

if length(sz) ~= 2
    error('Unsupported array size');
end

ain_is_logical = islogical(ain);
if ain_is_logical
    ain = double(ain);
end

aout = zeros(sz,class(ain));
ctr = mod(reshape(ctrpos, [1 2]) - 1, sz);

ctr_int = round(ctr);
ctr_dec = ctr - ctr_int;

if nargin > 2
    precision = varargin{1};
else
    precision = 1e-2;
end

use_dec = true;
if max(max(abs(ctr_dec))) < precision
    use_dec = false;
end

c2 = sz - ctr_int;


% First the integral shift
aout(1:c2(1),1:c2(2)) = ain(ctr_int(1)+1:end,ctr_int(2)+1:end);
aout(1:c2(1),c2(2)+1:end) = ain(ctr_int(1)+1:end,1:ctr_int(2));
aout(c2(1)+1:end,1:c2(2)) = ain(1:ctr_int(1),ctr_int(2)+1:end);
aout(c2(1)+1:end,c2(2)+1:end) = ain(1:ctr_int(1),1:ctr_int(2));

if use_dec
    faout = fftn(aout);
    [q1,q2] = ndgrid(-ceil(sz(1)/2):floor(sz(1)/2 - 1),-ceil(sz(2)/2):floor(sz(2)/2 - 1));
    q1 = fftshift(q1);
    q2 = fftshift(q2);
    aout = ifftn( exp(2i * pi * (q2 * ctr_dec(2) / sz(2) + q1 * ctr_dec(1) / sz(1))) .* faout);
end
    
if ain_is_logical
    aout = logical(aout);
end