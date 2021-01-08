%PLOT_OBJECTS plot reconstructed objects and layers  
% ** p              p structure
% ** use_display    if false, dont plot results on screen 
% *returns*
%  ++fig - image handle 


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


function [fig1, fig2] = plot_objects(p, use_display)

import math.*
import utils.*

% prepare quadratic phase for FP backpropagation and calculate the pixel
% size
[objpix, FP_pre_phase_factor] =  get_object_pixel_size(p);

% mask for backpropagation (FP only)
if p.fourier_ptycho && p.plot.filt
    for ii=1:p.numobjs
        ob_mask{ii} = ifftshift(filt2d_pad(p.object_size(ii,:), round(p.plot.FP_maskdim/p.dx_spec(1)*1.2), round(p.plot.FP_maskdim/p.dx_spec(1)), 'circ'));
    end
end

if length(unique(p.share_object_ID)) ~= length(p.object)
    % the GPU engine allows to modify the sharing within the engine and it
    % can cause inconsitencies during plotting
    utils.verbose(0,'Number of object does not correspond to the number of share_object_ID, resetting ... ')
    if length(p.object) == 1
        % assume shared scans
        p.share_object_ID(:) = 1; 
    elseif length(p.object) == length(p.share_object_ID)
        % assume unshared scans 
        p.share_object_ID(:) = 1:length(p.object); 
    else
        error('Correct settings of the shared objects could not be determined')
    end
end
   
count_plotobj = 1;

%modified by YJ for electron pty
if isfield(p,'beam_source') && strcmp(p.beam_source, 'electron')
    unitFactor = 0.1;
    unitLabel = 'nm';
else %X-ray 
    unitFactor = 1e6;
    unitLabel = '\mum';
end

for obmode = 1:p.object_modes
    for obnum = 1:p.numobjs
        
        % number of layers for multilayer object reconstruction 
        Nlayers = size(p.object{obnum},4); 
            
        % get the new object - 2D or 3D array
        ob_plot = p.object{obnum}(:,:,obmode,:); 
        
        % enforce update of the object_size
        object_size(obnum,:) = [size(ob_plot,1), size(ob_plot,2)]; 
        
        % get the actual reconstructed area and calculate the corresponding
        % mask
        ob_good_range = {p.asize(1)/2:object_size(obnum,1)-p.asize(1)/2, p.asize(2)/2:object_size(obnum,2)-p.asize(2)/2}; 
        plot_mask = false(object_size(obnum,:));
        plot_mask(ob_good_range{:},:) = true; 


        % remove phase offset and phase ramp (if requested)
        ob_plot = utils.stabilize_phase(ob_plot, 'weight', plot_mask, ...
            'remove_ramp', p.plot.remove_phase_ramp);

        
        if p.fourier_ptycho
            % propagate from lens plane to object plane
            ob_plot = ifft2(ifftshift(ob_plot.*ob_mask{obnum}))*p.object_size(obnum,1).*ifftshift(FP_pre_phase_factor{obnum});
        else
            
            if p.plot.show_layers
                if ~p.plot.show_layers_stack
                    % plot multiple layers next to each other
                    ob_plot = reshape(ob_plot,object_size(obnum,1), object_size(obnum,2)* Nlayers);
                    object_size(obnum,:) = size(ob_plot);
                    plot_mask = repmat(plot_mask,1,Nlayers);
                else
                    % 3D object for imagesc3D
                    ob_plot = squeeze(ob_plot);
                end

            else
                % plot single eDOF image 
                ob_plot = prod(ob_plot,4);  % show extended depth of focus images
                
            end
        end
        
        % apply apodization
        if p.plot.obj_apod
            try
                filt_size = [size(ob_good_range{1},2) size(ob_good_range{2},2)];
                ob_plot_size = [size(ob_plot,1),size(ob_plot,2)]; 
                ob_plot = ob_plot.*fftshift(utils.filt2d_pad(ob_plot_size, max(1,filt_size), max(1,filt_size-min(floor(filt_size.*0.05)))));
            catch
                utils.verbose(2, 'Failed to apply apodization.')
            end
        end
        
        % propagate object
        if p.plot.prop_obj ~= 0 
            ob_plot = utils.prop_free_nf(ob_plot, p.lambda, p.plot.prop_obj, objpix);
        end
        
        % get complex conjugate
        if p.plot.conjugate
            ob_plot = conj(ob_plot);
        end


        % precalculate absorption and phase
        absob = abs(ob_plot);
        phob = angle(ob_plot);

        
        %%%%%%%%%%%%%%%%%
        %%% AMPLITUDE %%%
        %%%%%%%%%%%%%%%%%

        % prepare figure handle for absorption images
        if ~use_display && count_plotobj == 1
            fig1 = plotting.smart_figure('Visible', 'off');
        else
            if count_plotobj == 1
                if p.plot.windowautopos && ~ishandle(1) % position it only if the window does not exist
                     fig1 = plotting.smart_figure(1);
                     set(gcf,'Outerposition',[ceil(p.plot.scrsz(4)/p.plot.horz_fact)+1 ceil(p.plot.scrsz(4)/2) ceil(p.plot.scrsz(4)/p.plot.horz_fact) ceil(p.plot.scrsz(4)/2)])    %[left, bottom, width, height
                 else
                     fig1 = plotting.smart_figure(1);
                 end
                clf;
            else
                set(groot,'CurrentFigure',fig1);
            end
        end

        ax_abs(count_plotobj)=subplot(p.plot.subplwinobj(1),p.plot.subplwinobj(2),count_plotobj);

        
        range = [min(p.positions([p.scanidxs{p.share_object_ID == obnum}],:)), ...
                 max(p.positions([p.scanidxs{p.share_object_ID == obnum}],:))];
        
        % FOV [xmin, ymin, xmax, ymax]
        good_fov(1) = -(p.object_size(obnum,1)/2 - p.asize(1)/2-range(1));
        good_fov(2) = -(p.object_size(obnum,2)/2 - p.asize(2)/2-range(2));
        good_fov(3) = good_fov(1) + range(3)-range(1);
        good_fov(4) = good_fov(2) + range(4)-range(2);
        good_fov = good_fov.*p.dx_spec([1,2,1,2])*unitFactor;
        fov_box = [good_fov(2),good_fov(1),good_fov(4)-good_fov(2), good_fov(3)-good_fov(1)]; % [xmin, ymin, W, H]  coordinates of the FOV box 
        
        % plot the absorption
        if ~p.plot.realaxes 
            plotting.imagesc3D(absob);
            if p.plot.fov_box && ~(p.plot.show_layers && ~p.plot.show_layers_stack && Nlayers > 1)
                rectangle('Position',[p.asize([2,1])/2 , p.object_size([2,1]) - p.asize([2,1])], 'EdgeColor', p.plot.fov_box_color)
            end
        else
            obj_ax = {([1 object_size(obnum,2)]-floor(object_size(obnum,2)/2)+1)*objpix(2)*unitFactor,([1 object_size(obnum,1)]-floor(object_size(obnum,1)/2)+1)*objpix(1)*unitFactor};
            plotting.imagesc3D(obj_ax{:},absob);
            xlabel(unitLabel)
            ylabel(unitLabel)
            if p.plot.fov_box && p.plot.show_layers && ~p.plot.show_layers_stack && Nlayers > 1
                % plot bar around each layer 
                for layer = 1:Nlayers
                    rectangle('Position', [good_fov(2) + (layer-(Nlayers+1)/2)*p.object_size(obnum,2).*p.dx_spec(1)*unitFactor ,good_fov(1), good_fov(4)-good_fov(2),good_fov(3)-good_fov(1)], 'EdgeColor', p.plot.fov_box_color)
                end
            elseif p.plot.fov_box
                rectangle('Position',fov_box, 'EdgeColor', p.plot.fov_box_color)
            end
        end
        
        % calculate a proper colorbar range
        try
            amp_range = sp_quantile(absob(plot_mask),[1e-4,1-1e-4], 10); 
        catch
            keyboard
        end
        if amp_range(1) < amp_range(2)
            caxis(amp_range);
        end
        
        colormap(bone(256)); colorbar
        axis image xy tight
        
        if check_option(p, 'show_only_FOV') && p.plot.realaxes
            axis([good_fov(2) good_fov(4) good_fov(1) good_fov(3)])
        elseif check_option(p, 'show_only_FOV') && ~p.plot.realaxes
            axis([p.asize(2)/2, object_size(obnum,2) - p.asize(2)/2, p.asize(1)/2, object_size(obnum,1) - p.asize(1)/2, ])
        end
        
        % prepare title strings
        if p.share_object
            title(sprintf('amplitude: %s %s', p.plot.obtitlestring, p.plot.extratitlestring),'interpreter','none');
        else
            title(sprintf('amplitude: %s %s', p.scan_str{obnum}, p.plot.extratitlestring),'interpreter','none');
        end
        
        
        %%%%%%%%%%%%%%%%%%%
        %%%%%% PHASE %%%%%%
        %%%%%%%%%%%%%%%%%%%
        
        % prepare figure handle for phase images
        if ~use_display && count_plotobj == 1
            fig2 = plotting.smart_figure('Visible', 'off');
        else
            if count_plotobj == 1
                 if p.plot.windowautopos && ~ishandle(2) % position it only if the window does not exist
                     fig2 = plotting.smart_figure(2);
                     set(gcf,'Outerposition',[1 ceil(p.plot.scrsz(4)/2) ceil(p.plot.scrsz(4)/p.plot.horz_fact) ceil(p.plot.scrsz(4)/2)])    %[left, bottom, width, height
                 else
                     fig2 = plotting.smart_figure(2);
                 end
                clf;
            else
                set(groot,'CurrentFigure',fig2);
            end
        end
        
        if check_option(p.plot, 'residua') && ~check_option(p, 'fourier_ptycho')
            % find residua to plot and avoid plotting residua in not illuminated regions 
            residues = plot_mask(2:end, 2:end) & (abs(utils.findresidues(ob_plot)) > 0.1); 
            [residues_ind{1}, residues_ind{2}] = find(residues); 
        end
        
        
        ax_phase(count_plotobj)=subplot(p.plot.subplwinobj(1),p.plot.subplwinobj(2),count_plotobj);
        
        % plot the phase
        if ~p.plot.realaxes
            plotting.imagesc3D(phob);
            if p.plot.fov_box && ~(p.plot.show_layers && ~p.plot.show_layers_stack && Nlayers > 1)
                rectangle('Position',[p.asize([2,1])/2 , p.object_size([2,1]) - p.asize([2,1])], 'EdgeColor', p.plot.fov_box_color)
            end
            if check_option(p.plot, 'residua')
                hold all
                plot(residues_ind{[2,1]},'or')
                hold off 
            end
        else
            plotting.imagesc3D(obj_ax{:},phob);
            xlabel(unitLabel)
            ylabel(unitLabel)
            if p.plot.fov_box && p.plot.show_layers && ~p.plot.show_layers_stack && Nlayers > 1
                % plot bar around each layer 
                for layer = 1:Nlayers
                    rectangle('Position', [good_fov(2) + (layer-(Nlayers+1)/2)*p.object_size(obnum,2).*p.dx_spec(1)*unitFactor ,good_fov(1), good_fov(4)-good_fov(2),good_fov(3)-good_fov(1)], 'EdgeColor', p.plot.fov_box_color)
                end
            elseif p.plot.fov_box
                rectangle('Position',fov_box, 'EdgeColor', p.plot.fov_box_color)
            end
            if check_option(p.plot, 'residua')
                hold all 
                plot((residues_ind{2}-size(ob_plot,2)/2)*objpix(2)*unitFactor,(residues_ind{1}-size(ob_plot,1)/2)*objpix(1)*unitFactor,'or')
                hold off 
            end
        end
        p_range = sp_quantile(phob(plot_mask),[1e-4,1-1e-4], 10); 
        if p_range(1) < p_range(2)
            caxis(p_range);
        end

        colormap(bone(256));colorbar 
        if p.share_object
            title(sprintf('phase: %s %s', p.plot.obtitlestring, p.plot.extratitlestring),'interpreter','none');
        else
            title(sprintf('phase: %s %s', p.scan_str{obnum}, p.plot.extratitlestring),'interpreter','none');
        end
        axis image xy tight
        if check_option(p, 'show_only_FOV') && p.plot.realaxes
            axis([-good_fov(2) good_fov(2) -good_fov(1) good_fov(1)])
        elseif check_option(p, 'show_only_FOV') && ~p.plot.realaxes
            axis([p.asize(2)/2, object_size(obnum,2) - p.asize(2)/2, p.asize(1)/2, object_size(obnum,1) - p.asize(1)/2, ])
        end
        count_plotobj = count_plotobj + 1;
    end
    
    if use_display 
        % link axes in case of zooming
        try linkaxes([ax_abs, ax_phase], 'xy'); end
    end
end

end



