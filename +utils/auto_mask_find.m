% mask = auto_mask_find(im,[<name>,<value>])
%
% im    Input complex valued image
%
% Optional parameters:
%
% margins   Two element array that indicates the (y,x) margins to exclude
%           from the edge of the mask window. For example to exclude the
%           noise around ptychography reconstructions, default 0.
% smoothing     Size of averaging window on the phase derivative, default
%               10.
% gradientrange Size of the histogram windown when selecting valid gradient
%               regions, in radians per pixel, default 1;
% show_bivariate    Show the bivariate histogram of the gradient, useful
%                   for debugging. Set to the number of figure you'd like
%                   it to appear.
%
%     Morphological operations to remove point details in the mask
%
% close_size    Size of closing window, removes dark bubbles from the mask, 
%               default 15.  ( = 1 for no effect)
% open_size     Size of opening window, removes bright bubbles from mask, 
%               default 120. ( = 1 for no effect)

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

function mask = auto_mask_find(im,varargin)
import plotting.franzmap

% Defaults
margin = [0 0];
smoothing = 10;
gradientrange = 1;
close_size = 15;
open_size = 120;
show_bivariate = 0;
zero_columns = [];

% parse the variable input arguments not handled by auto_mask_find
vararg = cell(0,0);                    
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch lower(name)
        case 'margin'
            margin = value;
        case 'smoothing'
            smoothing = value;
        case 'gradientrange'
            gradientrange = value;
        case 'close_size'
            close_size = value;
        case 'open_size'
            open_size = value;
        case 'show_bivariate'
            show_bivariate = value;    
        case 'zero_columns'
            zero_columns = value;
        otherwise
            vararg{end+1} = name;
            vararg{end+1} = value;
    end
end

% Some checks
if numel(margin)~= 2
    error('Margin variable should have two elements')
end

if close_size < 1
    error('erode_size must be an integer 1 or greater')
end

if open_size < 1
    error('erode_size must be an integer 1 or greater')
end

if gradientrange < 0
    error('gradientrange must be positive')
end

mask = true(size(im));


mask(1:1+margin(1),:) = false;
mask(end-margin(1):end,:) = false;
mask(:,1:1+margin(2)) = false;
mask(:,end-margin(2):end) = false;
if ~isempty(zero_columns)
    mask(:,zero_columns) = false;
end

% Compute phase gradient based on phasor
ph = exp(1i*angle(im));
[gx, gy] = gradient(ph);
gx = -real(1i*gx./ph);
gy = -real(1i*gy./ph);

kernel = ones(smoothing);
gx = conv2(gx,kernel,'same');
gy = conv2(gy,kernel,'same');

gaux(:,1) = gy(mask(:));
gaux(:,2) = gx(mask(:));

[N,C] = hist3_own(gaux,[100 100]);
[ny nx] = find(N == max(N(:)),1);

% masky = (gy>C{1}(ny-gradientrange))&(gy<C{1}(ny+gradientrange));
% maskx = (gx>C{2}(nx-gradientrange))&(gx<C{2}(nx+gradientrange));
masky = (gy>C{1}(ny)-gradientrange)&(gy<C{1}(ny)+gradientrange);
maskx = (gx>C{2}(nx)-gradientrange)&(gx<C{2}(nx)+gradientrange);

maskxy = maskx&masky;
% figure(1000); imagesc(masky); axis xy; colormap franzmap
% % Erosion
% erodemask = ones(erode_size);
% maskxy = erode_own(maskxy,erodemask);
% 
% % Dilation
% dilatemask = ones(dilate_size);
% maskxy = dilate_own(maskxy,dilatemask);

maskxy = close_own(maskxy,ones(close_size));
maskxy = open_own(maskxy,ones(open_size));

mask = mask&maskxy;

if show_bivariate > 0
    figure(show_bivariate); 
    imagesc(log10(N)); 
    colormap franzmap
end
end

function imout = erode_own(im,erodemask)
% My own erosion to avoid using Image Processing Toolbox
% Receives a binary image and kernel and performs erosion of the image
erodemask = erodemask/sum(erodemask(:));
imout = conv2(double(im),erodemask,'same');
imout = imout>0.99999;
end

function imout = dilate_own(im,dilatemask)
% My own dilation to avoid using Image Processing Toolbox
%  Receives a binary image and kernel and performs erosion of the image
dilatemask = dilatemask/sum(dilatemask(:));
imout = conv2(double(im),dilatemask,'same');
imout = imout>0;
end

function imout = open_own(im,openmask)
imout = dilate_own(erode_own(im,openmask),openmask);
end

function imout = close_own(im,closemask)
imout = erode_own(dilate_own(im,closemask),closemask);
end

function [histout, C] = hist3_own(gaux,bins)
eps = 0.001; % esther
min_gaux1 = min(gaux(:,1));
max_gaux1 = max(gaux(:,1));
inter_1   = (max_gaux1-min_gaux1)/bins(1);

min_gaux2 = min(gaux(:,2));
max_gaux2 = max(gaux(:,2));
inter_2   = (max_gaux2-min_gaux2)/bins(2);

indarray1 = floor( (1-eps)*bins(1)*( gaux(:,1)-min_gaux1 )./( max_gaux1-min_gaux1 ) + 1 );
indarray2 = floor( (1-eps)*bins(2)*( gaux(:,2)-min_gaux2 )./( max_gaux2-min_gaux2 ) + 1 );

histout = zeros(bins);

for ii = 1:numel(indarray1)
    histout(indarray1(ii),indarray2(ii)) = histout(indarray1(ii),indarray2(ii)) + 1;
end

C{1} = linspace(min_gaux1+inter_1/2,max_gaux1-inter_1/2,bins(1));
C{2} = linspace(min_gaux2+inter_2/2,max_gaux2-inter_2/2,bins(2));
end


