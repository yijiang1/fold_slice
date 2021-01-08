%PLOT_OBJECT_SPECTRUM 
% plot fourier transformation of the reconstructed object 
%
% ** p          p structure
% ** use_display    if false, dont plot results on screen 
%
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



function fig6 = plot_object_spectrum(p, use_display)

    count_plotobj = 1;

    [objpix] =  get_object_pixel_size(p);
    
    for obnum = 1:p.numobjs
        for obmode = 1:p.object_modes
                       
            ob_plot = p.object{obnum}(:,:,obmode,:);
            ob_plot = prod(ob_plot,4);  % make one eDoF image 
            if ~p.fourier_ptycho
                ob_plot = fftshift(fft2(ob_plot));
            end
            
            
            absob = abs(ob_plot);
            
            if ~use_display && count_plotobj == 1
                fig6 = plotting.smart_figure('Visible', 'off');
            else
                if count_plotobj == 1
                     if p.plot.windowautopos && ~ishandle(6) % position it only if the window does not exist
                         fig6 = plotting.smart_figure(6);
                         set(gcf,'Outerposition',[ceil(p.plot.scrsz(4)*2/p.plot.horz_fact)+1 ceil(p.plot.scrsz(4)/2) ceil(p.plot.scrsz(4)/p.plot.horz_fact) ceil(p.plot.scrsz(4)/2)])    %[left, bottom, width, height
                     else
                         fig6 = plotting.smart_figure(6);
                     end
                else
                    set(groot,'CurrentFigure',fig6);
                end
            end
            
            ax_abs(count_plotobj)=subplot(p.plot.subplwinobj(1),p.plot.subplwinobj(2),count_plotobj);
            
            good_fov(1) = (p.object_size(obnum,1)/2 - p.asize(1)/2) .*p.dx_spec(1)*1e6;
            good_fov(2) = (p.object_size(obnum,2)/2 - p.asize(2)/2) .*p.dx_spec(2)*1e6;
            
            if ~p.plot.realaxes
                imagesc(log10(absob));
            else
                obj_ax = {([1 p.object_size(obnum,2)]-floor(p.object_size(obnum,2)/2)+1)*objpix(2)*1e6,([1 p.object_size(obnum,1)]-floor(p.object_size(obnum,1)/2)+1)*objpix(1)*1e6};
                imagesc(obj_ax{:},log10(absob));
                xlabel('\mum')
                ylabel('\mum')
            end

            colorbar
            axis image xy tight

            if p.share_object
                title(sprintf('log10 object spectrum: %s %s', p.plot.obtitlestring, p.plot.extratitlestring),'interpreter','none');
            else
                title(sprintf('log10 object spectrum: %s %s', p.scan_str{obnum}, p.plot.extratitlestring),'interpreter','none');
            end
           
            
            count_plotobj = count_plotobj + 1;
        end
    end
    
end