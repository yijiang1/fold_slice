% FOURIER_RING_CORRELATION simplified but faster version of the FRC code
%
% [score, object] =  fourier_ring_correlation(object_1, object_2, varargin)
% 
% ** object_1    array reconstructed object 
% ** object_2    array reconstructed object from an independend scan 
% ** varargin    see code for more details 


%     
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
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% for LSQ-ML method 
% M. Odstrcil, A. Menzel, M.G. Sicairos,  Iterative least-squares solver for generalized maximum-likelihood ptychography, Optics Express, 2018
% for OPRP method 
%  M. Odstrcil, P. Baksh, S. A. Boden, R. Card, J. E. Chad, J. G. Frey, W. S. Brocklesby,  "Ptychographic coherent diffractive imaging with orthogonal probe relaxation." Optics express 24.8 (2016): 8360-8369
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
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



function [score, object] =  fourier_ring_correlation(object_1, object_2, varargin)

import engines.GPU.shared.*
import engines.GPU.GPU_wrapper.*
import math.*
import utils.*
import plotting.*

par = inputParser;
par.addParameter('px_scale',  1 , @isnumeric )
par.addParameter('auto_crop',  false, @islogical )
par.addParameter('plot_results',  true , @islogical )  % use white background 
par.addParameter('smoothing',  0 , @isnumeric )  % smooth over N pixels 
par.addParameter('Nrings',  20 , @isnumeric )  % smooth over N pixels 
par.addParameter('crop',  0 , @isnumeric )  % crop image by N pixels 
par.addParameter('flip_horizontal',  false , @islogical )  %  flip second image horizontally
par.addParameter('fft_phase_removal_guess',  false , @islogical )  %  flip second image horizontally
par.addParameter('weights',  {} , @iscell )  % cell array of weights 
par.addParameter('find_shift',  true, @islogical )  % cell array of weights 

par.parse(varargin{:})
r = par.Results;

if r.flip_horizontal
    object_2 = fliplr(object_2);
end

Npix = min(size(object_1), size(object_2));

object = {object_1, object_2}; 

if r.crop > 0
    for ii = 1:2
        object{ii} = crop_pad(object{ii}, Npix-r.crop);
    end
    if ~isempty(r.weights)
        for ii = 1:2
            r.weights{ii} = crop_pad(r.weights{ii}, Npix-r.crop);
        end
    end
end

Npix = min(size(object{1}), size(object{2}));

for ii = 1:2
    object{ii} = object{ii} / mean(abs(object{1}(:)) );

    W{ii} = tukeywin(Npix(1), 0.2) .* tukeywin(Npix(2),0.2)'; 
    if ~isempty(r.weights)
        W{ii} = W{ii} .* single(r.weights{ii}); 
    end
end

score.shift = [0,0];


if r.find_shift
    for kk = 1:4
        Npix = size(object{1});
        [X,Y] = meshgrid(-Npix(2)/2+1:Npix(2)/2,-Npix(1)/2+1:Npix(1)/2);
        object{1} = utils.stabilize_phase(object{1}, object{2}, 'fourier_guess', r.fft_phase_removal_guess);

        for ii = 1:2
            phasor{ii} = object{ii} ./ (abs(object{ii}) + 1e-3*mean(abs(object{ii}(:))));
            fobject{ii} = fft2(single(W{ii}.*(phasor{ii}-mean(phasor{ii}(:))))); 
        end


        % high pass filter 
        Wf = Garray(fftshift(exp(- 1./ ((X.^2+Y.^2)/(Npix(1)/50)^2)))); 

        [output] = utils.dftregistration( Wf.* fobject{1}, Wf.* fobject{2},100);
        object{2} = imshift_fft(object{2}, output(4), output(3));    


        ROI = { (1+max(0,ceil(output(3)))):(Npix(1)+min(0, floor(output(3)))) , ...
                (1+max(0,ceil(output(4)))):(Npix(2)+min(0, floor(output(4))))};

        object{1} = object{1}(ROI{:});
        object{2} = object{2}(ROI{:});
        for j = 1:2
            W{j} = W{j}(ROI{:});
        end
        
        verbose(3,'Image shifted by %g %g px', output([4,3]))
        score.shift = score.shift + Ggather([output(4), output(3)]);

    %     subplot(1,2,1)
    %     plotting.imagesc3D(object{1}); axis off image 
    %     subplot(1,2,2)
    %     plotting.imagesc3D(object{2}); axis off image 
    %     drawnow 

        if all(abs(output(3:4)) < 0.5)
            break
        end
    end
end

Npix = size(object{1});

[object{1}] = utils.stabilize_phase(object{1}, object{2}, abs(object{2}), 'binning', 4 , 'fourier_guess', r.fft_phase_removal_guess);
    
if r.flip_horizontal
    score.shift(1) = -score.shift(1);
end

for ii = 1:2
    object{ii} = object{ii} ./ mean(abs(object{ii}(:))); 
end

W = sqrt(W{1} .* W{2}); 
ROI_compare = get_ROI(W>0.1*max(W(:)));   % compare only the reliable ROIs

W = tukeywin(length(ROI_compare{1}),0.2) .* tukeywin(length(ROI_compare{2}),0.2)';

for ii = 1:2
    fobject{ii} = fft2(W.*object{ii}(ROI_compare{:})); 
end



for ii = 1:2
    fobject{ii} = fftshift(fobject{ii});
    fobject_norm{ii} =  abs(fobject{ii}).^2;
end
fcorr = fobject{1} .* conj(fobject{2}); 

binning = ceil(Npix/2  / r.Nrings);
fcorr = conv2(fcorr, ones(binning) / prod(binning), 'same');
fcorr = fcorr(1:binning(1):end, 1:binning(2):end);
for ii = 1:2
    fobject_norm{ii} = conv2(fobject_norm{ii}, ones(binning) / prod(binning), 'same');
    fobject_norm{ii} = fobject_norm{ii}(1:binning(1):end, 1:binning(2):end);
end

Npix= size(fcorr);


x = single(-Npix(2)/2+0.5:Npix(2)/2-0.5)/(Npix(2)/2);
y = single(-Npix(1)/2+0.5:Npix(1)/2-0.5)/(Npix(1)/2);
if  length(r.px_scale) > 1 && r.px_scale(1) > r.px_scale(2)
    y = y .* r.px_scale(2) / r.px_scale(1); 
elseif length(r.px_scale) > 1 && r.px_scale(1) < r.px_scale(2)
    x = x .* r.px_scale(1) / r.px_scale(2); 
end
[X,Y] = meshgrid(x, y);
    

R_mat = sqrt(X.^2 + Y.^2);

Rmax = 0.98; 
R0 = 0.01;
R_all = linspace(R0, Rmax, min(r.Nrings, min(Npix)));

for ii = 1:(length(R_all)-1)
    R = R_all(ii);
   
    ring = find((R_mat > R_all(ii)) & (R_mat < R_all(ii+1)));

    fcorr_values{1}(ii) = abs(sum(fcorr(ring)) ./ sqrt(sum(fobject_norm{1}(ring)) .* sum(fobject_norm{2}(ring)))); 

    n_values(ii) =  length(ring); %  sum(ring(:));
  
end

spatial_freq = R_all(1:end-1);

n_values = n_values .* prod(binning);
% 1-bit curve 
T = (0.5+2.41./sqrt(n_values)) ./ (1.5+1.41./sqrt(n_values));
% 1/2 bit curve 
% T = (0.21+1.91 ./sqrt(n_values)) ./ (1.21+0.91./sqrt(n_values));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fcorr_values{1} = Ggather(fcorr_values{1});

AUC = nanmean(fcorr_values{1}); % area undear curve criterion 
score.AUC = AUC; 
score.thresh = T;
score.FRC = fcorr_values{1};


score.spatial_freq = spatial_freq;
score.SSNR = 2 * score.FRC ./ (1-score.FRC);

score.SNR_avg = nansum(score.SSNR .* spatial_freq) / sum(spatial_freq);


Ts = smooth(T);
if r.smoothing > 0
    score.FRC = imgaussfilt(score.FRC,r.smoothing);
end
[x0,y0,iout,jout] = intersections(spatial_freq,score.FRC,spatial_freq, Ts,false);

if all(score.FRC  >= Ts')
    score.resolution = 1;
elseif all(score.FRC  <= Ts')
    score.resolution = 0;
elseif any(x0 > 0.1)
    score.resolution = min(x0(x0 > 0.1));
else
    score.resolution = min(x0);
end



if r.plot_results
    subplot(1,2,1)
    hold all
    b = plot(spatial_freq,score.FRC+randn*0.1,'LineWidth', 2);
    h = plot(spatial_freq, Ts, 'k--', 'LineWidth', 2);
    plot(x0, y0, 'o')
    xlabel('Spatial frequency / Nyquist')
    % ylabel('FRC')
    ylabel(sprintf('Fourier ring correlation, AUC=%3.3g', AUC))
    hold off
    ylim([0,1])
    xlim([0,1])
    % r = vline(resolution, '-k');
    legend([b, h], 'FRC', '1 bit threshold','Location','Best');
    grid on 

    subplot(1,2,2)
    hold all 
    plot(score.spatial_freq, score.SSNR);
    try; vline(score.resolution); end
    hline(1)
    set(gca, 'yscale', 'log')
    hold off 
    grid on 
    ylabel(sprintf('Spectral SNR, SNR_{avg}=%3.3g', score.SNR_avg))


%     width = 10;
%     aspect_ratio=4/3;
%     height = width / aspect_ratio;
%     % set size of the resulting image 
%     set(gcf, 'PaperPosition', [1.5 1.5 width height]);
    plotting.suptitle(sprintf('Resolution=%.3gnm  AuC=%.3g', min(r.px_scale) / score.resolution * 1e9,AUC))

end


verbose(3,'AUC %g', score.AUC)
verbose(3,'SNR %g', score.SNR_avg)
verbose(3,'resolution %g (%g nm)', score.resolution, min(r.px_scale) / score.resolution*1e9)

try
    verbose('SSNR 0.1 %g  0.5  %g, 0.9 %g \n', log10(quantile(score.SSNR, [0.1, 0.5, 0.9] )))
end 


    


end



