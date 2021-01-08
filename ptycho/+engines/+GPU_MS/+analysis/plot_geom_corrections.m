% PLOT_GEOM_CORRECTIONS plot position refinement statistics - position errors, directions and weights 
%
% plot_geom_corrections(self, mode, object, iter, par, cache)
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** mode      structure containing reconstruction parameters related to the selected incoherent mode
% ** object    cell of arrays, reconstructed object 
% ** iter      current iteration number 
% ** par       structure containing parameters for the engines 
% ** cache     structure with precalculated values to avoid unnecessary overhead
%

% FUNCTION plot_geom_corrections(self, mode, object, iter, par, cache)
% plot positiones updates 

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


function plot_geom_corrections(self, mode, object, iter, par, cache)
    
    import engines.GPU_MS.GPU_wrapper.*
    import math.*
    import plotting.*
    import utils.*
    
    pos = mode.probe_positions;     
    pos_0 = self.probe_positions_0;
    

    Nplots = 4*(iter >= par.probe_position_search && ~isempty(par.probe_geometry_model)) ...
        + (iter >= par.estimate_NF_distance) + ...
          (iter >= par.detector_rotation_search) + ...
          (iter >= par.detector_scale_search);
    

    if ~ishandle(16165)
        plotting.smart_figure(16165)
        set(gcf,'Outerposition',[100 100 Nplots*330 400])    %[left, bottom, width, height
    else
        plotting.smart_figure(16165)
    end
  
    plot_id = 0; 

    if iter >= par.probe_position_search && ~isempty(par.probe_geometry_model)
        clf()
        subplot(1,Nplots,1)
        hold all
        plot(mode.scales , '-'); axis tight 
		ylabel('Relative pixel scaling correction [-]')
		xlabel('Iteration')
        hold off 
        grid on 
        title('Scales')
        subplot(1,Nplots,2)
        plot(mode.rotation , '-'); axis tight 
		ylabel('Rotation  [deg]')
        title('Rotation')
		xlabel('Iteration')
        grid on 
        subplot(1,Nplots,3)
        plot(mode.shear , '-'); axis tight 
		ylabel('Shear  [deg]')
        title('Shear')
		xlabel('Iteration')
        grid on 
        subplot(1,Nplots,4)
        plot(mode.asymmetry*100 , '-'); axis tight 
		ylabel('Asymmetry [%]')
        title('Asymmetry')
		xlabel('Iteration')
        grid on 
        plot_id = 4; 
    end
    
    if iter >= par.estimate_NF_distance
        subplot(1,Nplots,plot_id+1)
        plot(mode.distances * 1e6 , '-'); 
        axis tight
        grid on 
		ylabel('Propagation distance [um]')
        title('Nearfield propagation distance')
		xlabel('Iteration')
        plot_id = plot_id + 1; 
    end
    
    if iter >= par.detector_rotation_search
        subplot(1,Nplots,plot_id+1)
        plot(mode.probe_rotation,'-'); 
        axis tight
        grid on 
		ylabel('Detector rotation angle [deg]')
        title('Detector rotation')
		xlabel('Iteration')
        plot_id = plot_id + 1; 
    end
    
    if iter >= par.detector_scale_search
        subplot(1,Nplots,plot_id+1)
        plot((1+mode.probe_scale_upd),'-'); 
        axis tight
        grid on 
		ylabel('Detector optimal scaling [-]')
        title('Relative pixel scale')
		xlabel('Iteration')
        plot_id = plot_id + 1; 
    end
    
    plotting.suptitle('Evolution of geometry parameters')

    %modified by YJ: remove check_option(par, 'probe_geometry_model') to
    %plot position correction even without geom refinement
    %if iter >= par.probe_position_search && check_option(par, 'probe_geometry_model')
	if iter >= par.probe_position_search
        %modified by YJ for electron pty
        if isfield(par,'beam_source') && strcmp(par.beam_source, 'electron')
            unitFactor = 1;
            scaleFactor = 0.1;
            unitLabel = 'A';
        else %X-ray 
            unitFactor = 1e9;
            scaleFactor = 1e6;
            unitLabel = 'nm';
        end
        % substract the geometry model to show only residuum
        pos_err = pos -  mode.probe_positions_model ; 

        % subtract average error per scan 
        for kk = 1:par.Nscans
           ind = self.reconstruct_ind{kk};
           pos_err(ind,:)  =  pos_err(ind,:) - mean(pos_err(ind,:)); 
        end        
        pos = pos+ self.Np_o([2,1])/2;

        marker_colors = {'r', 'b', 'g', 'k'};
        scale = self.pixel_size*scaleFactor; 
        plotting.smart_figure(455454)
        clf()
        subplot(2,2,1)
        aobject = angle(object);
        range = sp_quantile(aobject(cache.object_ROI{:}), [1e-3, 1-1e-3],10);
        aobject = (aobject  - range(1)) / (range(2) - range(1));
        grids = {(-ceil(self.Np_o(2)/2):ceil(self.Np_o(2)/2)-1)*scale(2), ...
                 (-ceil(self.Np_o(1)/2):ceil(self.Np_o(1)/2)-1)*scale(1)}; 
        imagesc(grids{:}, aobject, [-2, 1]); % reduce contrast 
        colormap bone 
        axis xy
        hold on 
        if isfield(par,'beam_source') && strcmp(par.beam_source, 'electron')
            ylabel('Position [nm]')
        else
            ylabel('Position [\mum]')
        end
        pos_scales = (pos-self.Np_o([2,1])/2) .* scale([2,1]); 
        for i = 1:length(self.reconstruct_ind)
            id = self.reconstruct_ind{i};
            if any(mode.probe_positions_weight)
                % plot importance 
                scatter(pos_scales(id,1), pos_scales(id,2), max(mode.probe_positions_weight(id,:),[],2)*20, marker_colors{1+mod(i,4)})
            end
            mean_err = mean(std(pos_err)); 
            range = max(pos) - min(pos); 
            up = 0.02 * min(range) / mean_err; 
            rounding_order = 10^floor(log10(up)); 
            up = ceil(up / rounding_order)*rounding_order; 
            quiver( pos_scales(id,1),  pos_scales(id,2),  scale(1)*pos_err(id,1)*up,  scale(2)*pos_err(id,2)*up, 0, marker_colors{1+mod(i,4)})
        end
        hold off 
        axis equal xy tight 
        range =  [min(pos_scales(:,1)), max(pos_scales(:,1)), min(pos_scales(:,2)), max(pos_scales(:,2))];
        axis(range)
        title(sprintf('Position errors, upscaled %ix', up))

        subplot(2,2,3)
        plot(mean(mode.probe_positions_weight,2), 'b.-')
        ylim([0, max(mean(mode.probe_positions_weight,2))])
        hold all
        for i = 1:length(self.reconstruct_ind)
            vline(self.reconstruct_ind{i}(end), '-r')
        end
        hold off 
        axis tight 
        ylabel('Importance weights')
        xlabel('Position #')   
        title('Relative importance weights for geometry model')

        subplot(2,2,2)
        yyaxis left 
        h = plot(pos_err(:,1), 'w.');
        axis tight 
        ylabel('Position error [px]')
        yyaxis right 
        plot(pos_err(:,1)*self.pixel_size(2)*unitFactor, 'b.-')
        axis tight 
        xlabel('Position #')
        ylabel(strcat('Position error [',unitLabel,']'))

        hold all
        for i = 1:length(self.reconstruct_ind)
            vline(self.reconstruct_ind{i}(end), '-r')
        end
        hold off 
        title( 'Horizontal')
        grid on 
        %legend({sprintf('STD=%3.2g nm', std(pos_err(:,1)*self.pixel_size(2)*1e9) )})
        legend({sprintf(strcat('STD=%3.2g ',unitLabel), std(pos_err(:,1)*self.pixel_size(2)*unitFactor) )})

        subplot(2,2,4)
        yyaxis left 
        h = plot(pos_err(:,2), 'w.');
        ylabel('Position error [px]')
        axis tight 
        yyaxis right 
        plot(pos_err(:,2)*self.pixel_size(1)*unitFactor, 'b.-')
        axis tight
        %ylabel('Position error [nm]')
        ylabel(strcat('Position error [',unitLabel,']'))

        xlabel('Position #')
        hold all
        title( 'Vertical')
        legend({sprintf(strcat('STD=%3.2g ',unitLabel), std(pos_err(:,2)*self.pixel_size(1)*unitFactor) )})
        grid on 

        for i = 1:length(self.reconstruct_ind)
            vline(self.reconstruct_ind{i}(end), '-r')
        end
        hold off
        plotting.suptitle('Random position errors after subtraction of geometry model')

        try
            if length(self.reconstruct_ind) == 2 && verbose()  > 1 && length(self.reconstruct_ind{1}) == length(self.reconstruct_ind{2})
                disp('Correlation between two scans')
                corr( pos_err(self.reconstruct_ind{1},:), pos_err(self.reconstruct_ind{2},:) )
            end
        end
    end
end