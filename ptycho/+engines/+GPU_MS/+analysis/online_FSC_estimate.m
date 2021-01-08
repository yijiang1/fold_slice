% ONLINE_FSC_ESTIMATE online estimation of the fourier shell correlation curve to estimation of optimal convergence
%  compare two scans and estimate FSC and other statistics 
%
% score = online_FSC_estimate(self, par, cache, score_0, iter)
%
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** par       structure containing parameters for the engines 
% ** cache     structure with precalculated values to avoid unnecessary overhead
% ** score_0      [] or a structure with outputs from previous online estimation of FSC curve
%
% returns:
% ++ score      structure with outputs from online estimation of FSC curve


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


function score = online_FSC_estimate(self, par, cache, score_0, iter)
    import engines.GPU_MS.GPU_wrapper.*
    import math.*
    import utils.*
    import plotting.*
    import engines.GPU_MS.*

    if check_option(self, 'object_orig')
        self.object{end+1,1} = cat(3,self.object_orig{1,:}); 
    end
    
    compared_indices =  (2:size(self.object,1))-1; 
    
    % take product of the reconstructed images, eDOF
    %% refererene image 
    selected_ROI = cache.object_ROI; 
    selected_ROI{2} = selected_ROI{2}(ceil(end/10):floor(end*9/10));
    obj{1} = cat(3, self.object{1,:});
    obj{1} = Garray(obj{1}); 
    
    for ll = compared_indices
        %% compared image 
        % take product of the reconstructed images, eDOF
        obj_compared = cat(3,self.object{ll+1,:});
        if ~isempty(score_0) && ~isempty(score_0{ll})
           obj_compared = imshift_fft(obj_compared, score_0{ll}.shift);
        end
        obj{2} = Garray(obj_compared);
        
        %% get at least some empirical esitmation of reliability -> for selection of compared ROI
        for kk = 1:2
            ind = [1,min(ll+1, length(cache.illum_sum_0))]; 
            W{kk} = cache.illum_sum_0{ind(kk)}(selected_ROI{:}); 
            W{kk} = W{kk} > 0.5*mean(W{kk}); 
            % W{kk} = imfill(W{kk}, 'holes');
            if any(W{kk}(:)==0) 
                Npix = size(W{kk}); 
                downscale = 10; 
                W{kk} = real(utils.interpolateFT(W{kk}, ceil(Npix / downscale)));
                try;  W{kk} = Garray(imerode( Ggather(W{kk})>0.1, strel('disk', ceil(self.Np_p(1)/8/downscale))));  end
                W{kk}  = (utils.imgaussfilt3_conv(W{kk}, mean(self.Np_p)/8/downscale)); 
                W{kk} = max(0,real(utils.interpolateFT(W{kk},Npix)));
           end
        end
        clear obj_0 
        Wshared = sqrt(W{1}.*W{2}); 
        for kk = 1:2
            W{kk} = Wshared; 
        end
        
        if size(obj{1},3) > 1 ||size(obj{2},3) > 1 
            Nl_shifts = 4;
        else
            Nl_shifts = 1;
        end
        
        for kk = 1:Nl_shifts
            
            for ii = 1:2
                Nlayers = size(obj{ii},3); 
                horiz_shifts = linspace(-(kk-1), (kk-1), Nlayers)'; 
                shift = [horiz_shifts, zeros(Nlayers,1)];
                if kk > 1 && ii == 1
                    shift = shift - score{ll,kk-1}.shift; 
                end
                % apply different shift on each layer -> minic rotation
                obj_tmp{ii} = prod(imshift_fft(obj{ii}, shift),3); 
                obj_tmp{ii} = obj_tmp{ii}(selected_ROI{:});
            end
            
            [score{ll,kk},obj_out] = analysis.fourier_ring_correlation(obj_tmp{:},...
                'smoothing', 1, 'crop', ceil(self.Np_p / 4) , 'plot_results', false, 'px_scale', self.pixel_size, 'weights', W);
        

            if ~isempty(score_0) && ~isempty(score_0{ll})
                score{ll,kk}.shift = score{ll,kk}.shift + score_0{ll}.shift ;
            end

            if ll == compared_indices(end) && verbose > 2
                plotting.smart_figure(2121)
                img = angle(cat(3,obj_out{:})); 
                plotting.imagesc3D(img); axis off image xy ; 
                caxis(Ggather(math.sp_quantile(img, [0.01, 0.99],10)))
                title('Aligned frames used for FSC estimation')
                drawnow 
            end

            % fprintf('========== total object shift ====== %g %g\n', score{end}.shift)
            score{ll,kk}.iter = iter; 
            score{ll,kk}.positions = self.modes{1}.probe_positions;
            score{ll,kk}.positions_0 = self.modes{1}.probe_positions_0;
            %score{ll,kk}.intensity = self.modes{1}.weights;
            score{ll,kk}.probe_fourier_shift =  self.modes{1}.probe_fourier_shift;
        end

    end
        

    
    
    plotting.smart_figure(4554)
    clf
    subplot(1,2,1)
    linestyle = {'-','--',':'}; 
    hold all
    for kk = 1:Nl_shifts
        for ll = compared_indices
            if isempty(score{ll,kk}); continue; end
            b(ll) = plot(score{ll,kk}.spatial_freq,score{ll,kk}.FRC,linestyle{1+mod(ll-1,end)},'LineWidth', 2);
            legend_names{ll} = sprintf('FRC scans 1 vs %i', ll+1);
        end
        h = plot(score{ll,1}.spatial_freq, score{ll,1}.thresh, 'k--', 'LineWidth', 2);
    end
    xlabel('Spatial frequency / Nyquist')
    ylabel(sprintf('Fourier ring correlation, AUC=%3.3g', score{ll,1}.AUC))
    hold off
    ylim([0,1])
    xlim([0,1])
    
    legend([b, h], legend_names{:}, '1 bit threshold','Location','Best');
    grid on 
    
    
    subplot(1,2,2)
    hold all 
    for kk = 1:Nl_shifts
        for ll = compared_indices
            if isempty(score{ll,kk}); continue; end
            score{ll,kk}.SSNR(~isfinite(score{ll}.SSNR) | score{ll,kk}.SSNR <= 0)  = nan; 
            plot(score{ll,kk}.spatial_freq, score{ll,kk}.SSNR);
        end
    end
    hline(1)
    set(gca, 'yscale', 'log')
    hold off 
    grid on 
    ylabel(sprintf('Spectral SNR, SNR_{avg}=%3.3g', score{ll,1}.SNR_avg))
    plotting.suptitle(sprintf('Resolution %3.3gnm', mean(self.pixel_size) / score{ll,1}.resolution * 1e9))



end

