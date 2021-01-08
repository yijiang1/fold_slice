% Integrates frames from a loopscan
%
% Syntax:
% [int] = integrate_frames(base_path,scan_num,plotfigure,det_num,savedata,maskfilename,masktype)
% Needed parameters: base_path (e.g. '~/Data10/')
%                    scan_num   (scan number)
% Optional parameters: plotfigure (figure number for final plot, 0 for no plotting, by default is 0)
%                      det_num (detector number, default 1)
%                      savedata (=1 to save data in 'analysis/integrated_frames/', default 0)
%                      maskfilename (valid mask file name. If empty [], no mask used)
%                      masktype (type of valid mask file name: 'bin' for binary or 'Oliver')
% 14-11-2012

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

function [int] = integrate_frames(base_path,scan_num,plotfigure,det_num,savedata,maskfile,masktype)
import beamline.pilatus_valid_pixel_roi
import io.image_read
import io.spec_read
import utils.compile_x12sa_filename

if exist('det_num') == 0
   det_num=1;
end
if exist('savedata') == 0
   savedata=0;
end
if exist('maskfile') == 0
   maskfile=[];
end
if exist('plotfigure') == 0
   plotfigure = 0;
end
if exist('masktype') == 0
   plotfigure ='Oliver';
end
%maskfile=[];
if savedata
   savefolder=[base_path 'analysis/integrated_frames/'];
   if savefolder ~= 7
      mkdir(savefolder)
   end
end

path=sprintf('%spilatus_%01d',base_path,det_num);
filename0=compile_x12sa_filename(scan_num,0);
data0=image_read(filename0);
img0=data0.data;
dims=size(img0);
S=spec_read(base_path,'ScanNr',scan_num);
num=size(S.bpm4i,1);

if ~isempty(maskfile)
    switch lower(masktype)
        case 'oliver'
            load(maskfile)
            valid_mask = pilatus_valid_pixel_roi(valid_mask,'RoiSize',size(img0));
            mask = zeros(size(img0));
            mask(valid_mask.indices) = 1;
        case 'bin'
            load(maskfile)
        otherwise
            disp('masktype unknown')
    end
else
    mask = ones(size(img0));
end

int=img0*0;
stack=zeros(dims(1),dims(2),num);

for jj=1:num
    filename=compile_x12sa_filename(scan_num,jj-1);
    data1=image_read(filename);
    img=data1.data.*mask;
    int=int+img;
    stack(:,:,jj)=img;
end

if plotfigure ~= 0
  figure(plotfigure)
  figure_position=[187   295   817   650];
  set(gcf,'Position',figure_position);
  imagesc(log10(int)); axis xy equal tight; colorbar; colormap jet
  title(sprintf('integrated frames S%05d',scan_num))
end
if savedata
  savefilename=sprintf('%s/S%05d_%01d_integrated_frames',savefolder,scan_num,det_num);
  save([savefilename '.mat'],'int')
  print('-f2','-djpeg','-r300',[ savefilename '.jpg'] );
  print('-f2','-depsc','-r1200',[savefilename '.eps'] );
end    