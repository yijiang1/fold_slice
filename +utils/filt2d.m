%FILT2D creates a 2d filter, based on fract_hanning
%
%   outputdim...                Size of the output array.
%   unmodsize...                Size of the central array containing no modulation. 
%   shape (optional)...         'rect' (default) or 'circ'
%   filter_type (optional)...   'hann' (default) or 'hamm', chebishev (only for unmodsize=0)
%
%   example: 
%       filt1 = filt2d(256,100,'circ','hann');
%       filt2 = filt2d(256,100);
%       filt3 = filt2d([512 420], [200 312]);
%   
%
%
% Adapted from fract_hanning:
%
% fract_hanning(outputdim,unmodsize)
% out = Square array containing a fractional separable Hanning window with
% DC in upper left corner.
% outputdim = size of the output array
% unmodsize = Size of the central array containing no modulation.
% Creates a square hanning window if unmodsize = 0 (or ommited), otherwise the output array 
% will contain an array of ones in the center and cosine modulation on the
% edges, the array of ones will have DC in upper left corner.

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

% License for fract_hanning:
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

function out = filt2d(outputdim, unmodsize, varargin)

if nargin > 2
    shape = varargin{1};
else
    shape = 'rect';
end

if nargin > 3
    filt_type = varargin{2};
else
    filt_type = 'hann';
end
    
if nargin == 1
    unmodsize = 0;
end

if length(outputdim)<2
    outputdim = [outputdim outputdim];
elseif length(outputdim)>2
    error('3D filters are not supported.')
end

if any(outputdim < unmodsize)
    error('Output dimension must be smaller or equal to size of unmodulated window'),
end

if unmodsize<0
    unmodsize = 0;
    warning('Specified unmodsize<0, setting unmodsize = 0')
end

if length(unmodsize)<2
    unmodsize = [unmodsize unmodsize];
elseif length(unmodsize)>2
    error('3D filters are not supported.')
end

switch lower(shape)
    case 'rect'
        N1 = [0:outputdim(2)-1];
        N2 = [0:outputdim(1)-1];
        [Nc,Nr] = meshgrid(N1,N2);
    case 'circ'
        assert(length(unique(outputdim))==1 && length(unique(unmodsize))==1, 'Option "circ" is supported for square arrays only.')
        N = [0:outputdim(1)-1];
        Nsz = (outputdim(1)-1)/2;
        xx = linspace(-Nsz,Nsz,outputdim(1));
        [x,y] = meshgrid(xx,xx);
        r = sqrt(x.^2 + y.^2);
        out = zeros(outputdim(1), outputdim(1));
end

if unmodsize == 0
    switch lower(filt_type)
        case 'hann'
            switch lower(shape) 
                case 'rect'
                    out = (1+cos(2*pi*Nc/outputdim(1))).*(1+cos(2*pi*Nr/outputdim(2)))/4;
                case 'circ'
                    out1d = (1+cos(2*pi*N/outputdim(1)))/2;
                    out1d = fftshift(out1d);
                    out(r<=Nsz) = interp1(xx,out1d,r(r<=Nsz));
                    out = ifftshift(out);
                otherwise
                    error('Unknown shape %s for filter %s', shape, filt_type);
            end
        case 'hamm'
            switch lower(shape) 
                case 'rect'
                    out = (0.54+0.46*cos(2*pi*Nc/(outputdim(1)-1))).*(0.54+0.46*cos(2*pi*Nr/(outputdim(2)-1)));
                case 'circ'
                    out1d = (0.54+0.46*cos(2*pi*N/(outputdim(1)-1)));
                    out1d = fftshift(out1d);
                    out(r<=Nsz) = interp1(xx,out1d,r(r<=Nsz));
                    out = ifftshift(out);                    
                    
                otherwise
                    error('Unknown shape %s for filter %s', shape, filt_type);
            end
        case 'chebyshev'
            switch lower(shape)
                case 'rect'
                    beta = cosh(1/outputdim(1)*acosh(10^5));
                    w1 = cos(outputdim(1)*acos(beta.*cos(pi*Nc/outputdim(1))))/(cosh(1/outputdim(1)*acos(beta)));
                    w1_fft = abs(fft(w1,[],2));
                    beta = cosh(1/outputdim(2)*acosh(10^5));
                    w2 = cos(outputdim(2)*acos(beta.*cos(pi*Nr/outputdim(2))))/(cosh(1/outputdim(2)*acos(beta)));
                    w2_fft = abs(fft(w2,[],1));
                    out = w2_fft.*w1_fft;
                otherwise
                    error('Unknown shape %s for filter %s', shape, filt_type);
                    
            end
        otherwise
            error('Unknown filter %s', filt_type);
    end

else
    switch lower(filt_type)
        case 'hann'
            switch lower(shape)
                case 'rect'
                    % Columns modulation
                    out = (1+cos(2*pi*(Nc- floor((unmodsize(2)-1)/2) )/(outputdim(2)+1-unmodsize(2))))/2;
                    if floor((unmodsize(2)-1)/2)>0
                        out(:,1:floor((unmodsize(2)-1)/2)) = 1;
                    end
                    out(:,floor((unmodsize(2)-1)/2) + outputdim(2)+3-unmodsize(2):length(N1)) = 1;
                    % Row modulation
                    out2 = (1+cos(2*pi*(Nr- floor((unmodsize(1)-1)/2) )/(outputdim(1)+1-unmodsize(1))))/2;
                    if floor((unmodsize(1)-1)/2)>0
                        out2(1:floor((unmodsize(1)-1)/2),:) = 1;
                    end
                    out2(floor((unmodsize(1)-1)/2) + outputdim(1)+3-unmodsize(1):length(N2),:) = 1;
                    
                    out = out.*out2;
                case 'circ'
                    out1d = (1+cos(2*pi*(N- floor((unmodsize(1)-1)/2) )/(outputdim(1)+1-unmodsize(1))))/2;
                    if floor((unmodsize(1)-1)/2)>0
                        out1d(1:floor((unmodsize(1)-1)/2)) = 1;
                    end
                    out1d(floor((unmodsize(1)-1)/2) + outputdim(1)+3-unmodsize(1):length(N)) = 1;
                    out1d = fftshift(out1d);
                    out(r<=Nsz) = interp1(xx,out1d,r(r<=Nsz));
                    out = ifftshift(out);
  
                otherwise
                    error('Unknown shape %s for filter %s', shape, filt_type);
            end
            
        case 'hamm'
            switch lower(shape)
                case 'rect'
                    % Columns modulation
                    out = (0.54+0.46*cos(2*pi*(Nc-floor((unmodsize(2)-1)/2))/(outputdim(2)-unmodsize(2))));
                    if floor((unmodsize(2)-1)/2)>0
                        out(:,1:floor((unmodsize(2)-1)/2)) = 1;
                    end
%                     keyboard
                    out(:,floor((unmodsize(2)-1)/2) + outputdim(2)+3-unmodsize(2):length(N1)) = 1;
                    % Row modulation
                    out2 = (0.54+0.46*cos(2*pi*(Nr-floor((unmodsize(1)-1)/2))/(outputdim(1)-unmodsize(1))));
                    if floor((unmodsize(1)-1)/2)>0
                        out2(1:floor((unmodsize(1)-1)/2),:) = 1;
                    end
                    out2(floor((unmodsize(1)-1)/2) + outputdim(1)+3-unmodsize(1):length(N2),:) = 1;
                    
                    out = out.*out2;
                case 'circ'
                    out1d = (0.54+0.46*cos(2*pi*(N-floor((unmodsize(1)-1)/2))/(outputdim(1)-unmodsize(1))));
                    if floor((unmodsize(1)-1)/2)>0
                        out1d(1:floor((unmodsize(1)-1)/2)) = 1;
                    end
                    out1d(floor((unmodsize(1)-1)/2) + outputdim(1)+3-unmodsize(1):length(N)) = 1;
                    out1d = fftshift(out1d);
                    out(r<=Nsz) = interp1(xx,out1d,r(r<=Nsz));
                    out = ifftshift(out);
                otherwise
                    error('Unknown shape %s for filter %s', shape, filt_type);

            end
            
        otherwise
            error('Unknown filter %s', filt_type);
    end
end  

end

