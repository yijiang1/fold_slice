% [PSD, freq] = power_spectral_density(img, varargin)
% Computes the power spectral density of the provided 3D image.
% Can handle non-cube arrays but assumes the voxel is isotropic
%
% Inputs:
% img               input image (2D or 3D)
%
% Parameters: 
% thickring         Normally the pixels get assigned to the closest integer pixel ring in Fourier domain. 
%                   With thickring the thickness of the rings is increased by
%                   thickring, so each ring gets more pixels and more statistics
% auto_binning      apply binning if dimensions are significanlty different along each axis
% mask              bool array equal to false for ignored pixels of the fft space  
%
% Outputs: 
% PSD               PSD curve values
% freq              normalized spatial frequencies to 1
%
% Example of use: 
% img = randn(512,512,512); 
% utils.power_spectral_density(img, 'thickring', 3);



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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PSD, freq] = power_spectral_density(img, air, varargin)
import math.isint
import utils.*

%%%%%%%%%%%%%%%%%%%%% PROCESS PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parser = inputParser;
parser.addParameter('thickring',  3 , @isnumeric )  % thick ring in Fourier domain
parser.addParameter('auto_binning',  true , @islogical )  % bin FRC before calculating rings, it makes calculations faster 
parser.addParameter('max_rings',  200 , @isnumeric )  % maximal number of rings if autobinning is used 
parser.addParameter('mask',  true, @islogical  )     % bool array, equal to false for ignored pixels of the fft space  
parser.addParameter('windowautopos',  true, @islogical  )     % automatically position plotted window
parser.addParameter('figure_id',  101, @isint)     % call figure(figure_id)


parser.parse(varargin{:})
param = parser.Results;


disp('Calculating PSD');

% remove masked values from consideration (i.e. for laminography)
Fimg = abs(bsxfun(@times,fftn(img) , param.mask+eps)).^2;

[ny,nx,nz] = size(img);
nmin = min(size(img));

% avoid edge artefacts 
img = img .* tukeywin(size(img,1),0.5) .* tukeywin(size(img,2),0.5)' .* reshape(tukeywin(size(img,3),0.5),1,1,[]);


thickring = param.thickring; 

if param.auto_binning 
    % bin the correlation values to speed up the following calculations 
    % find optimal binning to make the volumes roughly cubic 
    bin = ceil(thickring/4) * floor(size(img)/ nmin); 
    % avoid too large number of rings 
    bin = max(bin, floor(nmin ./ param.max_rings)); 
    
    if any(bin > 1)
        fprintf('Autobinning %ix%ix%i \n', bin)
        thickring = ceil(thickring / min(bin));
        % fftshift and crop the arrays to make their size dividable by binning number 
        if ismatrix(img); bin(3) = 1; end
        % force the binning to be centered 
        subgrid = {fftshift(ceil(bin(1)/2):(floor(ny/bin(1))*bin(1)-floor(bin(1)/2)-1)), ...
                   fftshift(ceil(bin(2)/2):(floor(nx/bin(2))*bin(2)-floor(bin(2)/2)-1)), ...
                   fftshift(ceil(bin(3)/2):(floor(nz/bin(3))*bin(3)-floor(bin(3)/2)-1))}; 
        if ismatrix(img); subgrid(3) = [] ; end 
        % binning makes the shell / ring calculations much faster
        Fimg = ifftshift(utils.binning_3D(Fimg(subgrid{:}), bin));
    end
else
    bin = 1;
end


[ny,nx,nz] = size(Fimg);
nmax = max([nx ny nz]);
nmin = min(size(img));


% empirically tested that thickring should be >=3 along the smallest axis to avoid FRC undesampling 
thickring = max(thickring, ceil(nmax/nmin)); 

param.thickring  = thickring;

rnyquist = floor(nmax/2);
freq = [0:rnyquist];

x = ifftshift([-fix(nx/2):ceil(nx/2)-1])*floor(nmax/2)/floor(nx/2);
y = ifftshift([-fix(ny/2):ceil(ny/2)-1])*floor(nmax/2)/floor(ny/2);
if nz ~= 1
    z = ifftshift([-fix(nz/2):ceil(nz/2)-1])*floor(nmax/2)/floor(nz/2);
else
    z = 0;
end
[X,Y,Z] = meshgrid(single(x),single(y),single(z));
index = (sqrt(X.^2+Y.^2+Z.^2));

clear X Y Z


Nr = length(freq); 
for ii = 1:Nr
    r = freq(ii);
    progressbar(ii,Nr)
    % calculate always thickring, min ring thickness is given by the smallest axis
    ind = index>=r-thickring/2 & index<=r+thickring/2 ; 
    ind = find(ind);  % find seems to be faster then indexing 
    auxFimg = Fimg(ind); 
    C(ii)  = sum(auxFimg);
    n(ii) = numel(ind);  % Number of points
end

n = n*prod(bin);  % account for larger number of elements in the binned voxels

PSD = abs(C) ./ n;
freq = freq/freq(end); 

figure(param.figure_id)
hold all 
plot(freq, PSD)
hold off 
set(gca, 'yscale', 'log')
ylabel('Power spectral density')
xlabel('Spatial frequency/Nyquist')
grid on 

    
end



