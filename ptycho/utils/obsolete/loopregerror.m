% Evaluate registration error in an angle series. Returns the error between
% subsequent images. The last image is evaluated against the first but
% flipped in x. Recieves image FT with DC in (1,1), the image should have
% had the center in center of array.
% filt_stackFT  FT of stack of images, previously filtered if needed
% deltastack    Estimates of positions
% xmax          Vector with positions of horizontal edges of rectangular 
%               window for registration
% ymax          Same as xmax but vertical edges of window

% Academic License Agreement
%
% Source Code
%
% Introduction 
% •	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
%   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the cSAXS 
%   ptychography MATLAB package computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").
%
% Terms and Conditions of the LICENSE
% 1.	LICENSOR grants to LICENSEE a royalty-free, non-exclusive license to use the PROGRAM for academic, non-commercial purposes, upon the terms and conditions 
%       hereinafter set out and until termination of this license as set forth below.
% 2.	LICENSEE acknowledges that the PROGRAM is a research tool still in the development stage. The PROGRAM is provided without any related services, improvements 
%       or warranties from LICENSOR and that the LICENSE is entered into in order to enable others to utilize the PROGRAM in their academic activities. It is the 
%       LICENSEE’s responsibility to ensure its proper use and the correctness of the results.”
% 3.	THE PROGRAM IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR 
%       A PARTICULAR PURPOSE AND NONINFRINGEMENT OF ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. IN NO EVENT SHALL THE LICENSOR, THE AUTHORS OR THE COPYRIGHT 
%       HOLDERS BE LIABLE FOR ANY CLAIM, DIRECT, INDIRECT OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION WITH THE PROGRAM OR THE USE 
%       OF THE PROGRAM OR OTHER DEALINGS IN THE PROGRAM.
% 4.	LICENSEE agrees that it will use the PROGRAM and any modifications, improvements, or derivatives of PROGRAM that LICENSEE may create (collectively, 
%       "IMPROVEMENTS") solely for academic, non-commercial purposes and that any copy of PROGRAM or derivatives thereof shall be distributed only under the same 
%       license as PROGRAM. The terms "academic, non-commercial", as used in this Agreement, mean academic or other scholarly research which (a) is not undertaken for 
%       profit, or (b) is not intended to produce works, services, or data for commercial use, or (c) is neither conducted, nor funded, by a person or an entity engaged 
%       in the commercial use, application or exploitation of works similar to the PROGRAM.
% 5.	LICENSEE agrees that it shall make the following acknowledgement in any publication resulting from the use of the PROGRAM or any translation of the code into 
%       another computing language:
%       "Data processing was carried out using the cSAXS ptychography MATLAB package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
%       Scherrer Institut, Switzerland."
%
% Additionally, any publication using the package, or any translation of the code into another computing language should cite for difference map:
% P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379–382 (2008). 
%   (doi: 10.1126/science.1158573),
% for maximum likelihood:
% P. Thibault and M. Guizar-Sicairos, Maximum-likelihood refinement for coherent diffractive imaging, New J. Phys. 14, 063004 (2012). 
%   (doi: 10.1088/1367-2630/14/6/063004),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
%   (doi: 10.1364/OE.24.029089).
% 6.	Except for the above-mentioned acknowledgment, LICENSEE shall not use the PROGRAM title or the names or logos of LICENSOR, nor any adaptation thereof, nor the 
%       names of any of its employees or laboratories, in any advertising, promotional or sales material without prior written consent obtained from LICENSOR in each case.
% 7.	Ownership of all rights, including copyright in the PROGRAM and in any material associated therewith, shall at all times remain with LICENSOR, and LICENSEE 
%       agrees to preserve same. LICENSEE agrees not to use any portion of the PROGRAM or of any IMPROVEMENTS in any machine-readable form outside the PROGRAM, nor to 
%       make any copies except for its internal use, without prior written consent of LICENSOR. LICENSEE agrees to place the following copyright notice on any such copies: 
%       © All rights reserved. PAUL SCHERRER INSTITUT, Switzerland, Laboratory for Macromolecules and Bioimaging, 2017. 
% 8.	The LICENSE shall not be construed to confer any rights upon LICENSEE by implication or otherwise except as specifically set forth herein.
% 9.	DISCLAIMER: LICENSEE shall be aware that Phase Focus Limited of Sheffield, UK has an international portfolio of patents and pending applications which relate 
%       to ptychography and that the PROGRAM may be capable of being used in circumstances which may fall within the claims of one or more of the Phase Focus patents, 
%       in particular of patent with international application number PCT/GB2005/001464. The LICENSOR explicitly declares not to indemnify the users of the software 
%       in case Phase Focus or any other third party will open a legal action against the LICENSEE due to the use of the program.
% 10.	This Agreement shall be governed by the material laws of Switzerland and any dispute arising out of this Agreement or use of the PROGRAM shall be brought before 
%       the courts of Zürich, Switzerland. 

function [errorreg grad] = loopregerror(deltastack,filt_stackFT,xmask,ymask)
import utils.dftregistration

[nr, nc, nim] = size(filt_stackFT);
errorreg = 0;

% Compute shifted stack FT
for ii = 1:nim-1,
    filt_stackFT(:,:,ii) = shiftpp3(filt_stackFT(:,:,ii),deltastack(1,ii),deltastack(2,ii));
end


filt_stack = ifft2(filt_stackFT);         % Compute the stack
filt_stack = filt_stack(ymask,xmask,:); % Registration window
filt_stackp1 = filt_stack(:,:,2:end);   % Slicing variable, takes more memory
                                            % necessary for parallelizaiton
for ii = 1:nim-1 % Error of ii vs ii+1
    output = dftregistration(filt_stack(:,:,ii),filt_stackp1(:,:,ii),0);
    errorreg(ii) = output(1);
end

% Compute first vs last
output = dftregistration(filt_stack(:,:,end),fliplr(filt_stackp1(:,:,1)),0);
errorreg(nim) = output(1);

if nargout == 2,
    grad = deltastack*0;
    stepe = 1e1; % For forward derivative
    Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
    Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
    [Nc,Nr] = meshgrid(Nc,Nr);
    expx = exp(1i*2*pi*(stepe*Nc/nc));
    expy = exp(1i*2*pi*(stepe*Nr/nr));
    for ii = 2:nim-1,
        auxim = filt_stackFT(:,:,ii).*expx; % FT{fn(x-(xn+eps))}  
        auxim = ifft2(auxim);        % fn(x-(xn+eps)) Compute the image
        auxim = auxim(ymask,xmask,:);       % Registration window
        output = dftregistration(auxim,filt_stackp1(:,:,ii),0);
        errorregepsx = output(1);    % En(xn+eps)
        grad(2,ii) = grad(2,ii) + (errorregepsx - errorreg(ii))/stepe;
        
        output = dftregistration(auxim,filt_stack(:,:,ii-1),0);
        errorregepsx = output(1);    % En-1(xn+eps)
        grad(2,ii) = grad(2,ii) + (errorregepsx - errorreg(ii))/stepe; 
        
        auxim = filt_stackFT(:,:,ii).*expy; % FT{fn(y-(yn+eps))}  
        auxim = ifft2(auxim);        % fn(y-(yn+eps)) Compute the image
        auxim = auxim(ymask,xmask,:);       % Registration window
        output = dftregistration(auxim,filt_stackp1(:,:,ii),0);
        errorregepsy = output(1);    % En(yn+eps)
        grad(1,ii) = grad(1,ii) + (errorregepsy - errorreg(ii))/stepe;
        
        output = dftregistration(auxim,filt_stack(:,:,ii-1),0);
        errorregepsy = output(1);    % En-1(xn+eps)
        grad(1,ii) = grad(1,ii) + (errorregepsy - errorreg(ii))/stepe; 
    end
    % grad(1)
    clear ii,
    auxim = filt_stackFT(:,:,1).*expx; % FT{fn(x-(xn+eps))}  
    auxim = ifft2(auxim);        % fn(x-(xn+eps)) Compute the image
    auxim = auxim(ymask,xmask,:);       % Registration window
    output = dftregistration(auxim,filt_stackp1(:,:,1),0);
    errorregepsx = output(1);    % En(xn+eps)
    grad(2,1) = grad(2,1) + (errorregepsx - errorreg(1))/stepe;
     
    output = dftregistration(auxim,fliplr(filt_stack(:,:,nim)),0);
    errorregepsx = output(1);    % En-1(xn+eps)
    grad(2,1) = grad(2,1) + (errorregepsx - errorreg(1))/stepe; 
        
    auxim = filt_stackFT(:,:,1).*expy; % FT{fn(y-(yn+eps))}  
    auxim = ifft2(auxim);        % fn(y-(yn+eps)) Compute the image
    auxim = auxim(ymask,xmask,:);       % Registration window
    output = dftregistration(auxim,filt_stackp1(:,:,1),0);
    errorregepsy = output(1);    % En(yn+eps)
    grad(1,1) = grad(1,1) + (errorregepsy - errorreg(1))/stepe;
        
    output = dftregistration(fliplr(auxim),filt_stack(:,:,nim),0);
    errorregepsy = output(1);    % En-1(xn+eps)
    grad(1,1) = grad(1,1) + (errorregepsy - errorreg(1))/stepe; 
    
    %grad(nim)
    auxim = filt_stackFT(:,:,nim).*expx; % FT{fn(x-(xn+eps))}  
    auxim = ifft2(auxim);        % fn(x-(xn+eps)) Compute the image
    auxim = auxim(ymask,xmask,:);       % Registration window
    output = dftregistration(auxim,fliplr(filt_stackp1(:,:,1)),0);
    errorregepsx = output(1);    % En(xn+eps)
    grad(2,nim) = grad(2,nim) + (errorregepsx - errorreg(nim))/stepe;
        
    output = dftregistration(auxim,filt_stack(:,:,nim-1),0);
    errorregepsx = output(1);    % En-1(xn+eps)
    grad(2,nim) = grad(2,nim) + (errorregepsx - errorreg(nim))/stepe; 
      
    auxim = filt_stackFT(:,:,nim).*expy; % FT{fn(y-(yn+eps))}  
    auxim = ifft2(auxim);        % fn(y-(yn+eps)) Compute the image
    auxim = auxim(ymask,xmask,:);       % Registration window
    output = dftregistration(auxim,fliplr(filt_stackp1(:,:,1)),0);
    errorregepsy = output(1);    % En(yn+eps)
    grad(1,nim) = grad(1,nim) + (errorregepsy - errorreg(nim))/stepe;
        
    output = dftregistration(auxim,filt_stack(:,:,nim-1),0);
    errorregepsy = output(1);    % En-1(xn+eps)
    grad(1,nim) = grad(1,nim) + (errorregepsy - errorreg(nim))/stepe; 
end

errorreg = sum(errorreg),
    
%     
% for ii = 1:nim-1,
%         aux1 = ifft2(shiftpp3(filt_stackFT(:,:,ii),deltastack(1,ii),deltastack(2,ii)));
%         aux2 = ifft2(shiftpp3(filt_stackFT(:,:,ii+1),deltastack(1,ii+1),deltastack(2,ii+1)));
%         aux1 = aux1(ymask,xmask);
%         aux2 = aux2(ymask,xmask);
%         output = dftregistration(aux1,aux2,0);
%         errorreg = errorreg + output(1);
% end
% INCLUDE THE LAST TO FIRST ERROR