% This function (1) removes ramp and (2) unwrap

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

function object_phase_unwrap = fun_ramp_unwrap(object, asize)
import utils.auto_mask_find
import utils.goldsteinunwrap2
import utils.remove_linear_phase_smart
 
deltax = asize(1)/2;                     % From edge of region to edge of image in x
rx = [1+deltax : size(object,2)-deltax];
ry  = [asize(1)/2-50 : size(object,1)-asize(1)/2+50];        % Range in y

%% Remove ramp
%%%%%%%%%% Tweak automask parameters %%%%%%%%%%
smoothing = 25;     % Size of smoothing window
gradientrange = 4;  % Range of the phase gradient taken. In units of histogram bins
close_size = 15;    % Size of closing window, removes dark bubbles from the mask ( = 1 for no effect)
open_size = 110;    % Size of opening window, removes bright bubbles from mask ( = 1 for no effect)
sidex = deltax*4;
maskzero_columnrange = [1+sidex : size(object,2)-sidex];  % Columns to be forced so the mask is zero. e.g. = [200:500], could be useful for middle of cylinders

show_bivariate = 0;     % Show the gradient bivariate histogram, input figure here
flag_plot = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pixel = [deltax+5,ry(end)-5];

mask1 = auto_mask_find(object,'margin',asize/2,'smoothing',smoothing,...
'gradientrange',gradientrange,'close_size',close_size,'open_size',open_size,'show_bivariate',show_bivariate,'zero_columns',maskzero_columnrange);
close;

[object ph_err] = remove_linear_phase_smart(object,'mask',mask1);
object_phase = angle(object);
object_size = size(object);
 
if flag_plot
    figure(2017); clf
    
    subplot(2,2,1);
    imagesc(angle(object).*mask1);
    axis xy equal tight
    colormap bone
    title('auto mask');
    
    subplot(2,2,2);
    absob = angle(object);
    % imagesc(([1 object_size(2)]-floor(object_size(2)/2)+1)*dx_spec(2)*1e6,([1 object_size(1)]-floor(object_size(1)/2)+1)*dx_spec(1)*1e6,absob);
    imagesc(absob);
    axis xy equal tight
    colormap bone

    subplot(2,2,[3 4]);
    line_mid = floor(size(object,1)/2);
    plot(angle(object(line_mid,:))); grid on;
    title(sprintf('Line %d profile',line_mid))

    xlabel('\mum')
    ylabel('\mum')

    set(gca,'fontsize',10,'fontweight','bold');
end
        
%% Unwrap
sel = object_phase(ry,rx);
%sel_unwrap=unwrap(sel,[],2);  % Matlab's routine, fast but problematic
%with noise
sel_unwrap=goldsteinunwrap2(sel,0);  % Goldstein's method, very robust (gets slow if there are residues)
img_unwrap=object_phase;
img_unwrap(ry,rx) = sel_unwrap;
img_unwrap(ry,rx) = img_unwrap(ry,rx) - 2*pi*round(img_unwrap(pixel(2),pixel(1))/(2*pi));
object_phase_unwrap = img_unwrap;


fprintf('-- Done removing ramp and unwrapping.\n');


