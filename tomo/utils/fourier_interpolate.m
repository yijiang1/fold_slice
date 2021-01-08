
% function fourier_interpolate(val,theta,numord)
% This function uses a Fourier transform approximation to obtain the zero
% and first coefficient of a periodic function (a sine wave fit) if only data
% of half a period is available. Data within 0 and 180 degrees (inclusive)
% is given. Completes and takes fft of data, theta is only used for some
% checks, should be given in radians.
%
% a*sin(x*pi/180+b)+c
% Output phase is in radians
%
% It is meant to work with equally spaced angles in the range of 180
% degrees, including 180 degrees. It allows not exactly equally spaced but then
% accuracy cannot be guaranteed.
%
% Optional parameter numord can be used in order to retrieve higher orders
% of the fit


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


function fit = fourier_interpolate(val,theta,numord,sym_flip)
% sym_flip  For a function that is only measured from 0 to 180 then the
%           data is flipped mirrored in order to complete the period of 
%           the sinewave

% theta = [0:45:180]*pi/180;
% val = 5 + sin(theta+0.5);
if nargin<3
    numord = 1;
end

if ~exist('sym_flip')
    sym_flip = true;
end
 
% if any(theta)>1.05*pi;
%     error('This routine cannot account for angles larger than 180 at the moment')
% end

if ~all(size(theta)==size(val))
    error(['theta and val dont have the same number of elements'])
end
val = val(:);
theta = theta(:);
clear auxval
if sym_flip
    auxval = [val(:);-fliplr(val(2:end-1))+val(1)+val(end)];
else
    auxval = val;
end
auxfft = fft(auxval)/numel(auxval);
fit.a = abs(auxfft(2:1+numord)*2*1i);
fit.b = angle(auxfft(2:1+numord)*2*1i);
fit.c = auxfft(1);


% thetafine = linspace(0,pi,100);
% figure(1);
% plot(theta,val,'.-')
% hold on
% plot(thetafine,fit.a*sin(thetafine+fit.b)+fit.c,'-r')
% hold off
% 
% figure(2);
% plot(auxval,'.-');
