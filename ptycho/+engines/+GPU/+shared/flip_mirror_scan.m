%  self = flip_mirror_scan(self, align_objects )
% join two scans obtained at 0 and 180 degrees to get a better estiamte
% of geometry errors 
% 
% Inputs: 
%   self     main data structure 
%   align_objects     find optimal shifts between the scans to match the 0 and 180 deg scan
    
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



function self = flip_mirror_scan(self, align_objects)

    
    import engines.GPU.*
    import utils.*
    import engines.GPU.GPU_wrapper.*

    if nargin < 2
        align_objects = true;   % find optimal shifts between the scans to match the 0 and 180 deg scan 
    end
    
    verbose(1,'=== Flipping 2nd scan ==== ')
    
    assert(length(self.reconstruct_ind) == 2, 'Two mirrored (0vs180deg) scans are required')
    
    for layer = 1:size(self.object,2)
        self.object{2,layer} = fliplr(self.object{min(end,2),layer}); 
    end


    ind = self.reconstruct_ind{2};
    
    
    %% FLIP THE 2ND SCAN POSITIONS  
    % keep the position distance from the reconstruction edge 
    

    self.probe_positions_0 = flip_positions(self, self.probe_positions_0, ind);
    if ~isempty(self.probe_positions)
        self.probe_positions = flip_positions(self, self.probe_positions, ind);
    end

   
    %% flip other inputs
    for ii = 1:length(self.probe)
       self.probe{ii}(:,:,2,:) =  fliplr(self.probe{ii}(:,:,min(end,2),:));
    end
    
    if ~isempty(self.diffraction)
        self.diffraction(:,:,ind) = fliplr(self.diffraction(:,:,ind));
    end
    if ~isempty(self.mask)
        if size(self.mask,3) == 1
            self.mask = repmat(self.mask, 1, 1, self.Npos);
        end
        self.mask(:,:,ind) = fliplr(self.mask(:,:,ind));
    end
        
    if isfield(self, 'affine_matrix')
        % flip nondiagonal terms of the affine matrix for account for the flipping 
        self.affine_matrix{2} = self.affine_matrix{2} .* [1,-1; -1,1]; 
    end
    
    
    
   if align_objects
       % using crosscorrelation find optimal shift between objects  
        self = shared.align_objects(self); 
    else
                
        % move the probe back to the center of the asize 
        for ii = 1:2
            [cx, cy] = math.center(abs(self.probe{1}(:,:,ii,1)).^2);
            self.probe{1}(:,:,ii,:) = utils.imshift_fft(self.probe{1}(:,:,ii,:), -cx, -cy);
            self.object{ii} = utils.imshift_fft(self.object{ii}, -cx, -cy);
        end
   end 
 
 
    
    
end

function pos_0 = flip_positions(self, pos_0, ind)
    % flip positions so that the flipped object does not move after
    % reconstruction 
    pos = pos_0(ind,1); 
    offset = ceil(self.Np_o(2)/2-self.Np_p(2)/2); 
    pos  = -(pos + offset); 
    left_offset = -max(pos); 
    right_offset = self.Np_o(2) -max(pos)-self.Np_p(2) ; 
    pos = pos + right_offset - left_offset;
    pos  = pos - offset; 
    pos_0(ind,1) = pos ; 

end