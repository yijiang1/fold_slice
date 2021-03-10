% [resolution FSC T freq n stat] = fourier_shell_corr_3D_2e(img1,img2,param, varargin)
% Computes the Fourier shell correlation between img1 and img2. It can also
% compute the threshold function T. Images can be complex-valued.
% Can handle non-cube arrays but assumes the voxel is isotropic
% Modified by YJ for electron ptychography
%
% Inputs:
%     **img1, img2        Compared images 
%     **param             Structure containing parameters
% *optional*: 
%     **dispfsc = 1;      Display results
%     **SNRt = 0.5        Power SNR for threshold, popular options:
%                         SNRt = 0.5;      1 bit threshold for average
%                         SNRt = 0.2071;   1/2 bit threshold for average
%     **thickring         Normally the pixels get assigned to the closest integer pixel ring in Fourier domain. 
%                         With thickring the thickness of the rings is increased by
%                         thickring, so each ring gets more pixels and more statistics
%     **auto_thickring    do not calculate overlaps if thickring > 1 is used
%     **st_title          optional extra title in the plot 
%     **freq_thr =0.05    mimimal freq value above which the resolution is detected  
%     **show_fourier_corr show 2D Fourier correlation 
%     **mask              bool array equal to false for ignored pixels of the fft space  
% 
% returns: 
%     ++resolution        [min, max] resolution estimated from FSC curve
%     ++FSC               FSC curve values
%     ++T                 Threshold values
%     ++freq              spatial frequencies 
%     ++stat              stat - structure containing other statistics such as
%                         SSNR, area under FSC curve. average SNR, ....

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%|                                                                       |
%*-----------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language this copyright should be retained and the authors 
%   and institution should be acknowledged in written form. Additionally 
%   you should cite the publication most relevant for the implementation 
%   of this code, namely
%   Vila-Comamala et al. "Characterization of high-resolution diffractive 
%   X-ray optics by ptychographic coherent diffractive imaging," Opt. 
%   Express 19, 21333-21344 (2011).
%   
%   Note however that the most relevant citation for the theoretical
%   foundation of the FSC criteria we use here is 
%   M. van Heela, and M. Schatzb, "Fourier shell correlation threshold
%   criteria," Journal of Structural Biology 151, 250-262 (2005).
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [resolution FSC T freq n stat] = fourier_shell_corr_3D_2e(img1,img2,param, varargin)
import math.isint
import utils.*

%%%%%%%%%%%%%%%%%%%%% PROCESS PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fsc_tic = tic;
if nargin < 3
    param = struct();
end

parser = inputParser;
parser.addParameter('dispfsc',  true , @islogical )
parser.addParameter('dispsnr',  false , @islogical )   % show also signal to noise ratio

parser.addParameter('SNRt', 0.5 , @isnumeric )% SNRt = 0.2071 for 1/2 bit threshold for average of 2 images
                                              % SNRt = 0.5 for 1 bit threshold for average of 2 images
parser.addParameter('thickring',  0 , @isnumeric )  % thick ring in Fourier domain
parser.addParameter('auto_binning',  false , @islogical )  % bin FRC before calculating rings, it makes calculations faster 
parser.addParameter('max_rings',  200 , @isnumeric )  % maximal number of rings if autobinning is used 
parser.addParameter('st_title',  '' , @isstring )  % optional extra title 
parser.addParameter('freq_thr',  0.05 , @isnumeric ) % mimimal freq value where resolution is detected  
parser.addParameter('show_2D_fourier_corr',  false , @islogical )  % instead of rings, show rather 2D distribution of the Fourier correlation
parser.addParameter('pixel_size',  []  )  % size of pixel in angstrom
parser.addParameter('mask',  [], @(x)(isnumeric(x) || islogical(x))  )     %  array, equal to 0 for ignored pixels of the fft space  and 1 for rest 
parser.addParameter('windowautopos',  true, @islogical  )     % automatically position plotted window
parser.addParameter('xlabel_type',  'nyquist', @(x)ismember(lower(x), {'nyquist', 'resolution'}))     % select X axis units 
parser.addParameter('figure_id',  100, @isint)     % call figure(figure_id)
parser.addParameter('clear_figure',  false, @islogical)     % clear figure before plotting 
parser.addParameter('out_fn',  [], @isstr)     % saving path for the image 
parser.addParameter('show_summary', true, @islogical) % show summary at the end

parser.parse(varargin{:})
r = parser.Results;

% load all to the param structure 
for name = fieldnames(r)'
    if ~isfield(param, name{1})  % prefer values in param structure 
        param.(name{1}) = r.(name{1});
    end
end

if isempty(param.pixel_size)
    warning('Pixel size not specified. Please use param.pixel_size. \n');
    param.pixel_size = nan; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create an example
% A = 3;
% img1 = rand(100,100,100);
% img2 = img1 + A*rand(100,100,100);
% img1 = img1 + A*rand(100,100,100);
% dispfsc = 1;
% SNRt = 1/A^2;

if any(size(img1) ~= size(img2))
    error('Images must be the same size')
end

    
[ny,nx,nz] = size(img1);
nmin = min(size(img1));


utils.verbose(2,'Calculating FSC');

% remove masked values from consideration (i.e. for laminography)
F1 = fftn(img1); 
F2 = fftn(img2); 
if ~isempty( param.mask)
    F1 = bsxfun(@times,F1 , param.mask+eps);
    F2 = bsxfun(@times,F2 , param.mask+eps);
end
F1cF2  = F1 .* conj(F2); 
F1 = abs(F1).^2; 
F2 = abs(F2).^2; 

[ny,nx,nz] = size(img1);
nmin = min(size(img1));

thickring = param.thickring; 

if param.auto_binning 
    % bin the correlation values to speed up the following calculations 
    % find optimal binning to make the volumes roughly cubic 
    bin = ceil(thickring/4) * floor(size(img1)/ nmin); 
    % avoid too large number of rings 
    bin = max(bin, floor(nmin ./ param.max_rings)); 
    
    if any(bin > 1)
        utils.verbose(1,'Autobinning %ix%ix%i', bin)
        thickring = ceil(thickring / min(bin));
        % fftshift and crop the arrays to make their size dividable by binning number 
        if ismatrix(img1); bin(3) = 1; end
        % force the binning to be centered 
        subgrid = {fftshift(ceil(bin(1)/2):(floor(ny/bin(1))*bin(1)-floor(bin(1)/2)-1)), ...
                   fftshift(ceil(bin(2)/2):(floor(nx/bin(2))*bin(2)-floor(bin(2)/2)-1)), ...
                   fftshift(ceil(bin(3)/2):(floor(nz/bin(3))*bin(3)-floor(bin(3)/2)-1))}; 
        if ismatrix(img1); subgrid(3) = [] ; end 
        % binning makes the shell / ring calculations much faster
        F1 = ifftshift(utils.binning_3D(F1(subgrid{:}), bin));
        F2 = ifftshift(utils.binning_3D(F2(subgrid{:}), bin));
        F1cF2 = ifftshift(utils.binning_3D(F1cF2(subgrid{:}), bin));
    end
else
    bin = 1;
end


[ny,nx,nz] = size(F1);
nmax = max([nx ny nz]);
nmin = min(size(img1));


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

% deal with asymmetric pixel size in case of 2D FRC
if  length(param.pixel_size) == 2
    if param.pixel_size(1) > param.pixel_size(2)
        y = y .* param.pixel_size(2) / param.pixel_size(1); 
    else
        x = x .* param.pixel_size(1) / param.pixel_size(2); 
    end
    param.pixel_size = min(param.pixel_size);  % FSC will be now calculated up to the maximal radius given by the smallest pixel size
end
                                
                                
[X,Y,Z] = meshgrid(single(x),single(y),single(z));
index = (sqrt(X.^2+Y.^2+Z.^2));

clear X Y Z


Nr = length(freq); 
for ii = 1:Nr
    r = freq(ii);
    if utils.verbose>2
        progressbar(ii,Nr)
    end
    % calculate always thickring, min ring thickness is given by the smallest axis
    ind = index>=r-thickring/2 & index<=r+thickring/2 ; 
    ind = find(ind);  % find seems to be faster then indexing 
    auxF1 = F1(ind); 
    auxF2 = F2(ind); 
    auxF1cF2 = F1cF2(ind); 
    C(ii)  = sum(auxF1cF2);
    C1(ii) = sum(auxF1);
    C2(ii) = sum(auxF2);
    n(ii) = numel(ind);  % Number of points
end


FSC = abs(C)./(sqrt(C1.*C2));
n = n*prod(bin);  % account for larger number of elements in the binned voxels

T = (  param.SNRt + 2*sqrt(param.SNRt)./sqrt(n+eps) + 1./sqrt(n)  )./...
    (  param.SNRt + 2*sqrt(param.SNRt)./sqrt(n+eps) + 1  );

freq_fine = 0:1e-3:max(freq);
freq_fine_normal = freq_fine/max(freq);

FSC_fine = max(0,interpn(freq, FSC, freq_fine, 'spline')); % spline, linear
T_fine = interpn(freq, T, freq_fine, 'spline');


idx_intersect = abs(FSC_fine-T_fine)<2e-4;
intersect_array = FSC_fine(idx_intersect);
range = freq_fine_normal(idx_intersect);
if length(range)<1
    range = [0 1];
    intersect_array = [1 1];
end


%%%%%% CALCULATE STATISTICS %%%%%%%%%%%%%%
pixel = param.pixel_size; % angstrom


range_start = range(find(range>param.freq_thr, 1, 'first'));
if isempty(range_start)
    range_start = range(1);
end

resolution = [pixel/range_start, pixel/range(end)];
fsc_mean = mean(FSC)/pixel;

% calculate SNR:  Huang, Xiaojing, et al. "Signal-to-noise and radiation exposure considerations in conventional and diffraction x-ray microscopy." Optics express 17.16 (2009): 13541-13553.
SSNR = 2 * FSC ./ (1-FSC);   % spectral signal to noise ratio 
SNR_avg = nansum(SSNR .* freq) / sum(freq);  % average SNR (should correspond to signal^2 / noise^2 )

st_title_full = sprintf('%s \n Pixel size %.3f A\n FSC: thickring %d, intersect (%.3f, %.3f) \n Resolution (%.3f, %.3f) A \n Area under FSC = %.3f, <FSC(A)> = %.3f SNR_avg=%.3f', ...
    param.st_title, pixel, param.thickring, range_start, range(end), pixel/range_start, pixel/range(end), mean(FSC), fsc_mean, SNR_avg);

if param.show_summary
    utils.verbose(1,['== FSC report: ==' st_title_full])
end

stat.fsc_mean = mean(FSC);
stat.fsc_mean_1A = fsc_mean;% Area in inverse angstrom
stat.SNR_avg = SNR_avg; 
stat.FSC = FSC;
stat.threshold = T;


%%%%% PLOT FOURIER SHELL CORRELATION %%%%%%%%%%%%%%%%%

if param.dispfsc
    fontsize = 12; % font size
    plotting.smart_figure(param.figure_id); 
    if param.dispsnr
        subplot(1,2,1); 
    end 
        
    if param.clear_figure; cla ; end 
    hold all
    plot(freq/freq(end), FSC, '-','linewidth',2);
    plot(freq/freq(end), T, 'r','linewidth',2);
    plot(range, intersect_array, 'go','markersize',6,'MarkerFaceColor','none','linewidth',2);
    grid on
    axis([0 1 0 1]);
    hold off 
    switch param.SNRt
        case 0.2071 , legend('FSC','1/2 bit threshold');
        case 0.5,     legend('FSC','1 bit threshold');
        otherwise ,   legend('FSC',['Threshold SNR = ' num2str(param.SNRt)]);
    end       
    
    set(gca,'fontweight','bold','fontsize',fontsize,'xtick',[0:0.1:1],'ytick',[0:0.1:1]);

    switch lower(param.xlabel_type)
        case 'nyquist'
            xlabel('Spatial frequency/Nyquist')
        case 'resolution'
            xaxis = [0:0.1:1];
            ticks = 1./xaxis* param.pixel_size;
            for i = 1:length(ticks)
                order = floor(log10(ticks(i)))-1;
                tick = round(ticks(i)/10^order)*10^order;
                if isnan(tick)
                    tick = [];
                end
                XTickLabel{i} = tick;
            end
            set(gca,'XTickLabel',XTickLabel);
            xlabel(gca, ['Half-period resolution [A]'])
    end
    if param.windowautopos
        win_size = [800 600]; 
        screensize = get( groot, 'Screensize' );
        set(gcf,'Outerposition',[100 min(270,screensize(4)-win_size(2)) win_size]);  %[left, bottom, width, height]
    end
    ylabel('Fourier shell correlation')
    title(st_title_full,'interpreter','none')
end



%%%%% PLOT SIGNAL TO NOISE RATIO %%%%%%%%%%%%%%%%%
%%  only approximation of SNR, dont use in publications 
if param.dispsnr
    fontsize = 12; % font size
    if param.dispfsc
        subplot(1,2,2); 
    else
        figure(param.figure_id); 
    end 
    if param.clear_figure; cla ; end 
    
    hold all
    plot(freq/freq(end), SSNR, '-','linewidth',2);
    grid on
    xlim([0 1]);
    set(gca, 'yscale', 'log')
    hold off 
     
    set(gca,'fontweight','bold','fontsize',fontsize,'xtick',[0:0.1:1]);

    switch lower(param.xlabel_type)
        case 'nyquist'
            xlabel('Spatial frequency/Nyquist')
        case 'resolution'
            xaxis = [0:0.1:1];
            ticks = 1./xaxis* param.pixel_size;
            for i = 1:length(ticks)
                order = floor(log10(ticks(i)))-1;
                tick = round(ticks(i)/10^order)*10^order;
                if isnan(tick)
                    tick = [];
                end
                XTickLabel{i} = tick;
            end
            set(gca,'XTickLabel',XTickLabel);
            xlabel(gca, ['Half-period resolution [A]'])
    end
    ylabel('Spectral signal to noise ratio')
    if ~(param.dispfsc)
        title(st_title_full,'interpreter','none')
    end
end

if ~isempty(param.out_fn)
    utils.verbose(1,'saving %s',param.out_fn);
    print('-djpeg','-r300',param.out_fn);
end

if utils.verbose>2
    toc(fsc_tic)
end

if param.show_2D_fourier_corr && nz == 1
    %% show 2D fourier correlation 
    C  = F1.*conj(F2);
    C1 = abs(F1).^2;
    C2 = abs(F2).^2;
    Nwin = 40;
    kernel = gausswin(Nwin, 3*nmax/ny) .* gausswin(Nwin, 3*nmax/nx)'; 
    C = conv2(fftshift(C),kernel,'same'); 
    C1 = conv2(fftshift(C1),kernel,'same'); 
    C2 = conv2(fftshift(C2),kernel,'same'); 

    figure(323)
    imagesc(abs(C) ./ sqrt(C1 .* C2))
    axis off square 
    colorbar 
    title('2D fourier correlation')
end
end



