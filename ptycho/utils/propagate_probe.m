% propagate_probe.m
% Warning: Currently working only for square pixels

%close all   % Recommended if improfile will be used (there is a bug with cursor positioning otherwise)

import utils.*

% Loading reconstruction
scan = 581;
fext = '.h5';
filename=['*_recons' fext];
file = dir(fullfile('~/Data10/analysis/', utils.compile_x12sa_dirname(scan), filename)); 
if length(file) > 1 ; warning('Multiple files are matching'); end
% load the last reconstruction in the provided scan number 
filename = fullfile(file(end).folder, file(end).name); 
title_str = file(end).name;
io.load_ptycho_recons(filename);

lambda=p.lambda;        % wavelength [m]
dis=p.z;                         % sample-detector distance [m]
asize=p.asize(1);
pixsize=p.dx_spec(1);     % pixel size of final reconstruction [m]

% Axial propagation parameters - Cut through (y,z)
rangez = [-0.7e-3 0.7e-3];      % Range in meters, currently at zero
step_num = 300;             % Number of steps

% Propagate to one plane
prop_dis= 1.5e-3;   % chosen propagation distance

% Display option
disp = 'hsv';        % Either 'hsv', 'amp', 'phase'

x = [-asize/2 asize/2]*pixsize;
%%%
%%% Code starts
%%%
scrsz = get(0,'ScreenSize');
mask=ones(asize,asize) - (abs(probe)==0);
%probe_corr=rmphaseramp(probe,mask);
probe_corr=remove_linearphase_v2(probe,mask,20);
%probe=probe_corr;
probe_hsv(:,:,3)=abs(probe)/max(max(abs(probe)));
probe_hsv(:,:,2)=ones(asize,asize);
probe_hsv(:,:,1)=angle(probe)/(2*pi)+0.5;
figure(1); 
clf
switch disp
    case 'hsv'
        imagesc(x*1e6,x*1e6,hsv2rgb(probe_hsv)); 
    case 'amp'
        imagesc(x*1e6,x*1e6,abs(probe));
        colormap bone
    case 'phase'
        imagesc(x*1e6,x*1e6,angle(probe));
        colormap bone
end
axis xy equal tight; 
title(['Probe ' title_str], 'interpreter', 'none');
xlabel('x [\mum]')
ylabel('y [\mum]')
set(gcf,'Outerposition',[1 1 500 500])
%
propdists = linspace(rangez(1),rangez(2),step_num);

back_propag_all = prop_free_nf(probe, lambda, propdists, pixsize);
propag=squeeze(back_propag_all(asize/2,:,:)); 
max_int = max(max(abs(back_propag_all).^2)); 

  

propag_hsv=zeros(asize,step_num,3);
propag_hsv(:,:,3)=abs(propag)/max(max(abs(propag)));
propag_hsv(:,:,2)=ones(asize,step_num);
propag_hsv(:,:,1)=angle(propag)/(2*pi)+0.5;
figure(8); 
clf
switch disp
    case 'hsv'
        imagesc(propdists*1e3,x*1e6,hsv2rgb(propag_hsv));
    case 'amp'
        imagesc(propdists*1e3,x*1e6,abs(propag));
        colormap bone
    case 'phase'
        imagesc(propdists*1e3,x*1e6,angle(propag));
        colormap bone
end

title(['Axial propagation ' title_str], 'interpreter', 'none');
xlabel('z [mm]')
ylabel('x [\mum]')
set(gcf,'Outerposition',[1 scrsz(4)-480 1000 480])


back = prop_free_nf(probe, lambda, prop_dis, pixsize);
back_hsv=zeros(asize,asize,3);
back_hsv(:,:,3)=abs(back)/max(max(abs(back)));
back_hsv(:,:,2)=ones(asize,asize);
back_hsv(:,:,1)=angle(back)/(2*pi)+0.5;
figure(2); 
clf
switch disp
    case 'hsv'
        imagesc(x*1e6,x*1e6,hsv2rgb(back_hsv)); 
    case 'amp'
        imagesc(x*1e6,x*1e6,abs(back)); 
        colormap bone
    case 'phase'
        imagesc(x*1e6,x*1e6,angle(back)); 
        colormap bone
end
axis xy equal tight;
title(['Probe propagated to ' num2str(prop_dis*1e3) ' mm ' title_str],'interpreter', 'none');
xlabel('x [\mum]')
ylabel('y [\mum]')
set(gcf,'Outerposition',[501 1 500 500])

utils.verbose(0, '==============  Focus positions: %g um =============', 1e6*propdists(math.argmax(max_int)))

%probe=back;
% save('probe_from_S00108_128x128_test_0_recons_00_propag_3_mm.mat','probe');

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
