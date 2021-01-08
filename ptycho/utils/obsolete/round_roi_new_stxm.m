% Another script for STXM evaluation of round roi scans
% 
% Does not use Gaussians centered on grid points which leads to inaccuarate
% values at the center positions (due to adding up all the Gaussians)
% but rather uses the Matlab's "griddata" function.

clear all
import utils.compile_x12sa_filename
import utils.verbose

%% scan parameters
scannumber = 2370;

scans = ['S' num2str(scannumber,'%05d')];
asize = [128 128];
%pathdir = sprintf('/afs/psi.ch/project/cxs/CDI/cSAXS_sxdm_2009_05_tomo/analysis/%s/', scans);
pathdir = sprintf('../../analysis/%s/', scans);
datafile = [pathdir, sprintf('%s_data_%03dx%03d.mat', scans, asize(1), asize(2))];

scan_type = 'round_roi';
dr = 1.5e-6;
lx = 30e-6;
ly = 36e-6;
nth = 5;

dx_spec = [65.4e-9 65.4e-9];

export_figures_to_svg = 0;

%% create positions

rmax = sqrt((lx/2)^2 + (ly/2)^2);
nr = 1 + floor(rmax/dr);
positions = [];
for ir=1:nr+1
    rr = ir*dr;
    dth = 2*pi / (nth*ir);
    for ith=0:nth*ir-1
        th = ith*dth;
        x1 = rr * cos(th);
        x2 = rr * sin(th);
        if( abs(x1) > lx/2 || (abs(x2) >= ly/2) )
            continue
        end
        positions(end+1,:) = [x1/dx_spec(2) x2/dx_spec(2)]; %#ok<AGROW>
    end
end
numpts = size(positions,1);

%% load data
if ~exist(datafile,'file') 
      while ~exist(compile_x12sa_filename(scannumber,numpts-1),'file');
          disp(['Waiting for scan ' scans ' to finish.'])
          pause(10)
      end
     verbose(1); core.prepare_data_2d(asize, numpts, scannumber, [100,268],'','','',147576);
end

load(datafile)



%% intialize variables for STXM analysis and upsampling

upsample =4; 
ndx = ceil(upsample*(2+lx/dr));
ndy = ceil(upsample*(2+ly/dr));
[xx,yy] = meshgrid((1:ndx) - ndx/2, (1:ndy) - ndy/2);

pos = positions;
pos(:,1) = upsample*pos(:,1)*dx_spec(1)/dr;
pos(:,2) = upsample*pos(:,2)*dx_spec(2)/dr; 

trans = zeros(ndy,ndx);
dpcx = zeros(ndy,ndx);
dpcy = zeros(ndy,ndx);

% put the STXM results into linear vectors first
lin_trans = zeros(numpts,1);
lin_dpcx = zeros(numpts,1);
lin_dpcy = zeros(numpts,1);

%% STXM analysis loop

for ii=1:numpts
    s1 = sum(data(:,:,ii),1);
    s2 = sum(data(:,:,ii),2);
    Atrans  = sum(s1);
    Adpcx = sum(((1:asize(1))-asize(1)/2).*s1)/Atrans;
    Adpcy = sum(((1:asize(1))-asize(1)/2)'.*s2)/Atrans;
    lin_trans(ii) = Atrans;
    lin_dpcx(ii) = Adpcx;
    lin_dpcy(ii) = Adpcy;
end

%% upsample

trans=griddata(pos(:,1),pos(:,2),lin_trans,xx,yy,'linear');
dpcx_nearest=griddata(pos(:,1),pos(:,2),lin_dpcx,xx,yy,'nearest');
dpcx_linear =griddata(pos(:,1),pos(:,2),lin_dpcx,xx,yy,'linear');
%dpcy=griddata(pos(:,1),pos(:,2),lin_dpcy,xx,yy,'linear');

%I0 = mean(trans(10:20,3));
[n,x] = hist(trans(:),asize(1));
I0 = x(find(diff(n)>0,1,'last')+1);

 % Compute integrated linear attenuation coefficient
mu = -log(trans./ I0 );

%% plot

%figure(1)
% scatter(positions(:,1),positions(:,2),300,lin_dpcx,'filled');
% colormap bone(256)
 %colorbar
% title([scans ': scattered photons per diffraction pattern'])
% axis ij image off

% figure(2);
%imagesc(dpcx_nearest);
%axis image ij tight off
%colormap bone(256); colorbar
%title([scans ': scattered photons (upsampled by a factor of ' num2str(upsample) ', nearest neighbour interpolation)'])
%scalebar(dr,'fontsize',20, 'linewidth',10)

figure(1);
imagesc(dpcx_linear);
axis image ij tight off
colormap bone(256); colorbar
title([scans ': DPC x']) 
scalebar(upsample*dx_spec(1),'fontsize',20, 'linewidth',10)

figure(4)
imagesc(trans)
axis image ij tight off
colormap bone(256); colorbar
scalebar(upsample*dx_spec(1),'fontsize',20, 'linewidth',10)
title([scans ': transmission']) 

figure(5)
imagesc(mu)
axis image ij tight off
colormap bone(256); colorbar
scalebar(upsample*dx_spec(1),'fontsize',20, 'linewidth',10)

% 
% roi_half = floor(upsample*nr/sqrt(2));
% roi_mask = zeros(size(trans));
% roi_mask(size(trans,1)/2-roi_half+1:size(trans,1)/2+roi_half,...
%     size(trans,2)/2-roi_half+1:size(trans,2)/2+roi_half) =1;
% roi_mask(size(trans,1)/2-upsample+1:size(trans,1)/2+upsample,...
%     size(trans,2)/2-upsample+1:size(trans,2)/2+upsample) =0;
% 
% roi_ind = find(roi_mask == 1);

%% export figures

if export_figures_to_svg
saveas(1,fullfile(pathdir,'/S00327_dpcx_scatterplot.png'))
%plot2svg(fullfile(pathdir,'/S00327_dpcx_scatterplot.svg'),1)    
plot2svg(fullfile(pathdir,'/S00327_dpcx_interpol_nearest.svg'),2)
plot2svg(fullfile(pathdir,'/S00327_dpcx_interpol_linear.svg'),3)
plot2svg(fullfile(pathdir,'/S00327_trans_interpol_linear.svg'),4)
plot2svg(fullfile(pathdir,'/S00327_abs_mu_interpol_linear.svg'),5)
disp(['Image files saved.']);
end

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