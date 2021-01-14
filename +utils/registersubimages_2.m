function [subim1, subim2, delta, deltafine, regionsout] = registersubimages_2(img1, img2, x1, y1, x2, y2, upsamp, displ,Wfilt);
% Finds common subregions in img1 and img2 and computes registration (alignment)
% Does not perform global search, so an approximate initial estimate is
% necessary
%
% [subim1, subim2, delta, deltafine,regionsout] = registersubimages_2(img1, img2, x1, y1, x2, y2, upsamp, displ,Wfilt);
%
% Inputs
% img1 -    Reference image
% x1, y1 -  Vectors indicating the region of interest of img1, e.g. x1 =
%           [30:100]. If left empty the code starts with full window and 
%           then reduces to the common part
% img2, x2, y2 -    Same as 1 but for the image to align
% upsamp -  Fraction of a pixel for the fine search, default = 1
% displ     = 0  no information displayed (default)
%           = 1  to display text information
%           > 1  also shows images in figure(display)
% Wfilt     Fourier domain filter for registration (= 1 for no effect)
%
% Outputs 
% subim1, subim2 -  Common subimages taken from img1 and img2, subim2 is realigned subpixel
% delta          -  Registration values (dy,dx)
% deltafine      -  Only the subpixel components of the alignment
% regionsout     -  Structure with final regions for registration
%
% Please cite and acknowledge if you use it. Algorithm extended from that in 
% M. Guizar-Sicairos, S. T. Thurman and J. R. Fienup, "Efficient subpixel
% image registration algorithms," Opt. Lett. 33, 156 (2008).
%
% Report bugs to mguizar@gmail.com
%
% This code is intended for finding the same object in two different
% images. It is a local search, so the code will go with the first solution
% it finds, an approximately good initial estimate is required.
%
% This code uses integer pixel and subpixel cross-correlation to find the
% local solution. The intended purpose is for this code to be able to find
% two similar subimages in the same image so the global search would be
% inadequate. This code does not remove the global phase between the
% subimages.
%
% Keep your eyes open, this code has not been throughly debugged or tested
%
% History
%
% Allow for iterative refinement of chosen subwindow
% Allow for empty input for subwindow regions
% Manuel Guizar - January 22, 2013
%
% Modified on 24 Aug 2010 to return deltafine and receive a Fourier filter
% Also modified to have initial estimate with upsamp 2 and then round and
% so that the final registration is done with an upsam accuracy
%
% Modified on 25 Aug 2010. If coarse alignment exceeds image dimensions
% assume error in coarse alignment and perform only subpixel.
%
% Original code by Manuel Guizar - August 4, 2009

% Copyright (c) 2010, Manuel Guizar Sicairos, University of Rochester
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
import utils.dftregistration
import utils.shiftpp2

maxloops = 100;

if nargin<6,
    error('Not enough input arguments'),
elseif nargin == 6,
    upsamp = 1;
    displ = 0;
elseif nargin == 7,
    displ = 0;
elseif nargin > 9,
    error('Too many input arguments'),
end

if any(size(img1) ~= size(img2))
    error('Input images must be the same size')
end

[nc2,nr2] = size(img2);
if isempty(x1)
    x1 = [1:nr2];
end
if isempty(x2)
    x2 = x1;
end
if isempty(y1)
    y1 = [1:nc2];
end
if isempty(y2)
    y2 = y1;
end


if displ >= 1,
    display('Coarse alignment');
end
subim1 = img1(y1,x1);
S1 = fft2(subim1).*Wfilt;
delta = [0 0];
%%%%%%%%%%%%%%%%%%%%%%%%
%%% Single pixel alignment
%%%%%%%%%%%%%%%%%%%%%%%%
correction = 1;
counter = 0;
correprev = [0 0];

while correction == 1,
%     if (min(y2+delta(1))<1)||(min(x2+delta(2))<1)||(max(y2+delta(1))>nc2)||(max(x2+delta(2))>nr2)
%         warning('Coarse estimate exceeded array dimensions. This means that either the registration failed miserably or the register window does not exist in the second image.')
%         error('No registration was performed')
%         subim2 = img2(y2,x2); 
%         delta = [0 0]; 
%         deltafine = [0 0];
%         return
%     end
    %%% Check if window still exists in subimage 2 %%%%
    if min(x2+delta(2))<1
        if displ>0
            disp('Refining range in lower x')
        end
        x2 = x2(min(x2)-delta(2):end);
        x1 = x1(min(x1)-delta(2):end);
        subim1 = img1(y1,x1);
        S1 = fft2(subim1).*Wfilt;
    end
    if max(x2+delta(2))>nr2
        if displ>0
            disp('Refining range in upper x')
        end
        x2 = x2(1:max(x2)-delta(2));
        x1 = x1(1:max(x1)-delta(2));
        subim1 = img1(y1,x1);
        S1 = fft2(subim1).*Wfilt;
    end
    if min(y2+delta(1))<1
        if displ>0
            disp('Refining range in lower y')
        end
        y2 = y2(min(y2)-delta(1):end);
        y1 = y1(min(y1)-delta(1):end);
        subim1 = img1(y1,x1);
        S1 = fft2(subim1).*Wfilt;
    end
    if max(y2+delta(1))>nc2
        if displ>0
            disp('Refining range in upper y')
        end
        y2 = y2(1:max(y2)-delta(1));
        y1 = y1(1:max(y1)-delta(1));
        subim1 = img1(y1,x1);
        S1 = fft2(subim1).*Wfilt;
    end
    subim2 = img2(y2+delta(1),x2+delta(2)); 
    errors = dftregistration(S1,fft2(subim2),2);
    errors(3:4) = round(errors(3:4));
    
    if displ > 0
        display(['E = ' num2str(errors(1)) ', Correction displacement ('...
            num2str(errors(3)) ',' num2str(errors(4)) '), delta = '...
            num2str(delta) ])
    end
    if displ > 1
        figure(ceil(abs(displ)));
        subplot(1,2,1)
        imagesc((real(subim1)));
        title('Reference subimage')
        axis equal;
        axis xy;
        axis tight;
        colorbar;
        colormap bone
        subplot(1,2,2)
        imagesc((real(subim2)));
        title('Aligned subimage')
        axis equal;
        axis xy;
        axis tight;
        colorbar;
        colormap bone
    end
    if (abs(errors(3))==0)&&(abs(errors(4))==0)
        correction = 0;
    end
    delta = delta-errors(3:4);
    counter = counter + 1;
    if counter > maxloops
        if displ>0
            warning('Maximum number of iterations exceeded, revise the problem or increase number of iterations')
        end
        break;
    end
    
    if prod(([errors(3) errors(4)] == -correprev)*1.0) == 1,
        break; % Stuck going back and forth
    end
    correprev = [errors(3) errors(4)];
        
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Subpixel alignment
%%%%%%%%%%%%%%%%%%%%%%%%

%%% use shiftpp2 on the large image, then extract....   clever ah?
if upsamp > 1,
    if displ > 1,
        display('Fine alignment'),
    end
    correction = 1;
    deltafine = [0 0];
    damp = 1;
    counter = 0;
    while correction == 1;  
        subim2 = img2(y2+delta(1),x2+delta(2)); 
        errors = dftregistration(S1,fft2(subim2),upsamp);
        img2 = shiftpp2(img2,-damp*errors(3),-damp*errors(4)); %% Suboptimal, change to use a routine that receives FT data
        deltafine = deltafine+damp*errors(3:4);
        if displ > 0,
            display(['E = ' num2str(errors(1)) ', Refine correction ('...
            num2str(deltafine(1)) ',' num2str(deltafine(2)) ')' ])
        end 
        if displ > 1,
            figure(ceil(abs(displ)));
            subplot(1,2,1)
            imagesc((real(subim1)));
            title('Reference subimage')
            axis equal;
            axis xy;
            axis tight;
            colorbar;
            colormap bone
            subplot(1,2,2)
            imagesc((real(subim2)));
            title('Aligned subimage')
            axis equal;
            axis xy;
            axis tight;
            colorbar;
            colormap bone
        end
         
        if (abs(errors(3))<=1/upsamp)&&(abs(errors(4))<=1/upsamp)
            correction = 0;
        end
        counter = counter + 1;
        if counter > maxloops
            if displ>0
                warning('Maximum number of iterations exceeded in subpixel registration')
            end
            delta = delta-deltafine;
            regionsout.x1 = x1;
            regionsout.x2 = x2;
            regionsout.y1 = y1;
            regionsout.y2 = y2;
            return;
        end
    end
    delta = delta-deltafine; 
    regionsout.x1 = x1;
    regionsout.x2 = x2;
    regionsout.y1 = y1;
    regionsout.y2 = y2;
    %counter;
    if displ > 0
        disp(['Final registration values, delta = (' num2str(delta(1)) ',' num2str(delta(2)) ')'])
    end
end
    
    
    
    