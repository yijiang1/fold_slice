%  self = align_objects(self )
% align all provided object arrays and apply the same shift on the
% corresponding probe so that the scans can be used for shared object
% reconstructions 
% 
% Inputs: 
%   self     main data structure 
    
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



function self = align_objects(self)

    import engines.GPU.*
    import utils.*
    import engines.GPU.GPU_wrapper.*


    cache.skip_ind = [];
    [cache.oROI_s{1},cache.oROI{1}] = shared.find_reconstruction_ROI( self.probe_positions_0,self.Np_o, self.Np_p);
    cache.object_ROI = {ceil(self.Np_p(1)/2):self.Np_o(1)-ceil(self.Np_p(1)/2), ...
                        ceil(self.Np_p(2)/2):self.Np_o(2)-ceil(self.Np_p(2)/2)};

    Nobj = size(self.object,1); 
    shift = [0,0];   
    obj1 = prod(cat(3,self.object{1,:}),3);
    obj1 = Ggather(obj1(cache.object_ROI{:}));

    for ll = 1:Nobj
        obj2 = prod(cat(3,self.object{ll,:}),3);
        obj2 = Ggather(obj2(cache.object_ROI{:}));
        [score(ll), object_aligned{ll}] = analysis.fourier_ring_correlation(obj1,obj2,...
            'smoothing', 5, 'crop', self.Np_p/4, 'plot_results', false);    
        shift(ll,:) =  score(ll).shift; 
    end
    
    if any([score(:).AUC] < 0.3)
        warning('Alignment of scans %i vs scan %i propably failed, FRC resolution is %3.3g of Nyquist frequency limit', 1, ll, score.resolution)
        figure(56476)
        subplot(1,Nobj,1)
        plotting.imagesc3D(angle(object_aligned{1}{1})); axis off image
        grid on
        for ll = 1:Nobj
            subplot(1,Nobj,ll)
            plotting.imagesc3D(angle(object_aligned{ll}{2})); axis off image xy 
            grid on
        end
        plotting.suptitle('Objects after alignement, check visually the estimated alignement')
        colormap bone 
        disp('Estimated shifts between the aligned objects')
        disp(shift )
    end
    
    shift = shift - mean(shift,1);

    
    self.Np_o = self.Np_o + 2*ceil(max(abs(shift(:,[2,1])),[],1));
    for layer = 1:size(self.object,2)
        for ll = 1:Nobj
            % make them all the same size 
            self.object{ll,layer} = crop_pad(self.object{ll,layer}, self.Np_o, mean(self.object{ll,layer}(:))); 
            self.object{ll,layer} = utils.imshift_fft(self.object{ll,layer}, shift(ll,:)); 
        end
    end
    if utils.verbose > -1 && any(max(abs(shift) ./ self.Np_p([2,1])) > 0.25)
        warning off backtrace
        id = math.argmax(max(abs(shift),[],2));
        
        warning('Alignement of two mirrored scans resulted in maximal probe shift of %3.3g %3.3gpx, \n!! this is more than 50%% of the probe diameter !! \nconsider using better alignement between 0deg and 180deg scans ', -shift(id,1), shift(id,2))
        warning on backtrace
    end
    
    
    for ii = 1:length(self.probe)
        for ll = 1:Nobj
            % apply shift to the probe as well, if needed, replicate the
            % probe 
            probe{ii}(:,:,ll,:) = utils.imshift_fft(self.probe{ii}(:,:,min(end,ll),:), shift(ll,:)); 
        end
    end
    self.probe = probe ;

end