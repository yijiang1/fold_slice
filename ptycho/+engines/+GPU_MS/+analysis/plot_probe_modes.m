% PLOT_PROBE_MODES plot incoherent probe modes 
%
% plot_probe_modes(self, cache)
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** cache     structure with precalculated values to avoid unnecessary overhead
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
% for mixed coherent probe:
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


function plot_probe_modes(self, par)
    import engines.GPU_MS.GPU_wrapper.*
    import utils.*
    import math.*
    import plotting.*
    
    for i = 1:par.probe_modes
        power(i) = Ggather(mean2(abs(self.probe{i}(:,:,1)).^2));
    end
    power  = power / sum(power);

    grids = {(-ceil(self.Np_p(2)/2):ceil(self.Np_p(2)/2)-1)*self.pixel_size(2), ...
         (-ceil(self.Np_p(1)/2):ceil(self.Np_p(1)/2)-1)*self.pixel_size(1)}; 

    
    plotting.smart_figure(46456)
           
    for i = 1:par.probe_modes
        mode = Ggather(mean(mean(self.probe{i},3),4));
        
        ax(2*i-1)=subplot(2,par.probe_modes,i);
        RANGE = sp_quantile(abs(mode),[1e-3,1-5e-3], 4)';
        RANGE(2) = max(RANGE(2), RANGE(1)+1e-6);
        amode = abs(mode);
        imagesc3D(grids{:},amode)
        if diff(RANGE)>0;caxis(RANGE); end
        title(sprintf('Mode %i, P:%3.2g', i, power(i)))
        ylabel(sprintf('Amplitude - <%3.2g ; %3.2g>', RANGE))
        axis image xy
        colormap bone 
        set(gca,'TickLength',[0 0])     
        set(gca,'XTick',[],'YTick',[])
        
        
        ax(2*i)=subplot(2,par.probe_modes,par.probe_modes+i);
        arg = -angle(utils.stabilize_phase(mode));
        RANGE_arg = sp_quantile(arg,[1e-3,1-1e-3], 4)';
        imagesc3D(grids{:},arg )
        if diff(RANGE_arg)>0;  caxis((RANGE_arg')); end
        ylabel(sprintf('Phase - <%3.2g ; %3.2g>', RANGE_arg))
        axis image xy
        colormap bone 
        set(gca,'TickLength',[0 0])
        set(gca,'XTick',[],'YTick',[])
    end
    linkaxes(ax, 'xy')


    if par.probe_modes > par.Nscans  % dont run for multiscan 
       reconstruct_ind = [self.reconstruct_ind{:}];
        if par.variable_probe
            plotting.smart_figure(id+1)
            clf       
            power  = power / sum(power);
             for i = 1:length(self.probe)
                pos = self.probe{i}.probe_positions;
                pos = pos(:,[2,1]); 
                pos(:,1) = -pos(:,1);
                subplot(2,1,1)
                hold all
                plot(power(i))
                hold off 
                title('Variable incoherent probe')
                xlabel('Normalized mode power')
                subplot(2,par.probe_modes,par.probe_modes+i)
                scatter(pos(reconstruct_ind,:),gather(W(reconstruct_ind)), 20);
                axis off image 
            end
        end
        
    end

end
