% PLOT_RESULTS show current reconstruction and errors during ptychography
%
% plot_results(self, cache, par, fourier_error,probe_positions)
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** cache     structure with precalculated values to avoid unnecessary overhead
% ** par       structure containing parameters for the engines 
% ** fourier_error  array [Npos,1] containing evolution of reconstruction error 
% ** probe_positions  array [Npos,2] with probe positions for the main coherence mode 

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
% 


function plot_results(self, cache, par, fourier_error,probe_positions)


    import engines.GPU_MS.GPU_wrapper.*
    import utils.*
    import math.*
    import plotting.*

    likelihood = lower(par.likelihood);

    try    
        verbose(1,'Plotting ... ')
        Np_o = self.Np_o;
        Npos = length(probe_positions);  
        reconstruct_ind = [self.reconstruct_ind{:}];
        ind = find(any(~isnan(fourier_error),2));

        probe = Ggather(self.probe{1}(:,:,1,1)); 
        % show extended DoF projection through layers of first scan  
        object = prod(cat(3,self.object{1,:}),3);
        
        if par.fourier_ptycho
            object = fft2(fftshift(object)); 
        end
        object = object(cache.object_ROI{:});
     
        plotting.smart_figure(10)
        clf()
        number_iterations = size(fourier_error,1);
        ha = tight_subplot(2,2,[.01 .01],[.01 .01],[.01 .01]);
        axes(ha(1))

  
        pixel_size = self.pixel_size .* cosd(par.sample_rotation_angles([1,2]));
        params = {'scale', pixel_size,'enhance_contrast', true};
        probe_positions = probe_positions - repmat([mean(cache.object_ROI{2})-Np_o(2)/2, mean(cache.object_ROI{1})-Np_o(1)/2],Npos,1); 

        imagesc_hsv(object ,params{:});

        if ~par.fourier_ptycho
            % avoid plotting residua in not illuminated regions 
            resid_mask = cache.illum_sum_0{1}(cache.object_ROI{:})/ cache.MAX_ILLUM(1) > 0.1; 
            % find residua to plot 
            residues = resid_mask(2:end, 2:end) & (abs(utils.findresidues(object)) > 0.1); 

            [X,Y] = find(residues); 
            if length(probe_positions) < 2e3
                points = probe_positions(reconstruct_ind, :); 
                hold all 
                plot(points(:,1)*pixel_size(2)*1e6, points(:,2)*pixel_size(1)*1e6, '.w')
                plot((Y-size(object,2)/2)*pixel_size(2)*1e6,(X-size(object,1)/2)*pixel_size(1)*1e6,'ow')
                hold off 
            end
        end
        
        axis xy
        ylabel('Reconstruction in fake colors')
        axes(ha(3))
        probe = utils.prop_free_nf(probe, self.lambda, sum(self.z_distance(1:end-1))/2, self.pixel_size); 
        imagesc_hsv(probe, params{:} );
 
        
        axis xy
        ylabel('Contrast enhanced probe')
        subplot(2,2,2)
        fourier_error(fourier_error == 0) = nan;
        

        if strcmpi(likelihood, 'poisson')
            fourier_error = (bsxfun(@minus, fourier_error, fourier_error(1,:))); 
        end
         if ~isempty(ind)  %if there is somethign to plot 
            hold all
            plot(ind, fourier_error(ind,reconstruct_ind), '-')
            ind_missing = ~ismember(1:self.Npos, reconstruct_ind);
            if any(ind_missing)
                plot(ind, fourier_error(ind,ind_missing),  '--')
            end
            plot(ind, nanmean(fourier_error(ind,~ind_missing)'),'k', 'LineWidth', 3)
            plot(ind, nanmedian(fourier_error(ind,~ind_missing)'),'k--', 'LineWidth', 3)
            hold off 
            grid on 
            set(gca, 'xscale', 'log')
            if strcmpi(likelihood, 'L1')
                set(gca, 'yscale', 'log')
            end
            xlim([1, number_iterations])
            % ignore the first iteration error in plotting 
            try   ylim([min2(fourier_error(2:end,:)), max2(fourier_error(2:end,:))]); end

            switch likelihood
                case 'poisson', title('Relative neg-likelihood change');
                case  'l1', title('Fourier error'); 
            end
         end
        subplot(2,2,4)
        if length(ind) > 1
            err = fourier_error(ind(end) , reconstruct_ind) ;
            pos = pixel_size([2,1]).*probe_positions(reconstruct_ind,:);
            %% compatibility with the CPU code 
            pos(:,2) = -pos(:,2);
            show_spatial_distribution(Ggather(pos), Ggather(err), false, length(probe_positions) < 2e3)
            axis off equal
        end
        title('Spatial distribution of error')


    catch err
        warning('Error during plotting: %s', err.message)
        keyboard
        disp('plotting failed')
    end  

end
