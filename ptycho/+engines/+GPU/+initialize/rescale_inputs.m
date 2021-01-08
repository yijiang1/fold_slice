% RESCALE_INPUTS multigrid scheme method
%
% self = rescale_inputs(self, Np_p_new, rescale_data)
%
% 
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** Np_p_new  new size of the rescaled dataset 
% ** rescale_data  true -> rescale also data + mask + noise arrays 
%
% returns: 
% ++ self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..

    
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

function self = rescale_inputs(self, Np_p_new, rescale_data)

    import utils.*
    import math.*

    if isempty(Np_p_new)
        return 
    end
    
    scale = Np_p_new ./ self.Np_p; 
    
    if all(scale == 1)
        return
    end
    
    self.modes = [];

    if ~isempty(self.diffraction_deform_matrix)
        verbose(1,'Rescalling with diffraction_deform_matrix is not implemented')
        self.diffraction_deform_matrix = [];
    end
    
    self.Np_p = ceil(self.Np_p.*scale);
    self.Np_o = ceil(self.Np_o.*scale);

    if ~isempty(self.probe_positions) && any(self.probe_positions(:) ~= self.probe_positions_0(:))
        self.probe_positions = scale([2,1]) .* self.probe_positions;
    else
        self.probe_positions = []; 
    end
    self.probe_positions_0 = scale([2,1]) .* self.probe_positions_0;
    self.pixel_size = self.pixel_size ./ scale; 

    for i = 1:length(self.probe)
        % scale also intensity 
        self.probe{i} = interpolateFT_addnoise(self.probe{i}, self.Np_p, 1 )./prod(scale);
        %self.probe{i} = imresize(self.probe{i}, self.Np_p)./prod(scale);
    end

    for i = 1:numel(self.object)
        self.object{i} = interpolateFT_addnoise(self.object{i}, self.Np_o , 1);
        %self.object{i} = imresize(self.object{i}, self.Np_o);
    end
    
    if ~isempty(self.probe_support)  
        self.probe_support = interpolateFT( self.probe_support, self.Np_p) ;    
    end
    if ~isempty(self.probe_support_fft)  
        self.probe_support_fft = crop_pad( self.probe_support_fft, self.Np_p) ;    
    end
    
    if rescale_data
        self.diffraction = fftshift_2D(self.diffraction);
        self.diffraction = crop_pad(self.diffraction, self.Np_p);
        self.diffraction = ifftshift_2D(self.diffraction);
        
        if ~isempty(self.mask)
            fill_value = 0.9;  % fill value in case of ptychographic "super resolution"
            self.mask = fftshift_2D(self.mask);
            self.mask = crop_pad(self.mask, self.Np_p, fill_value);
            self.mask = ifftshift_2D(self.mask);
        end
        if ~isempty(self.noise)
            self.noise = fftshift_2D(self.noise);
            self.noise = crop_pad(self.noise, self.Np_p);
            self.noise = ifftshift_2D(self.noise);

        end
    end
    
        

end

function [ imout ] = interpolateFT_addnoise(im,outsize, noise_amplitude)
    % fill the empty regions in the FFT interpolated data by some weak
    % random noise 
    import math.fftshift_2D
    import math.ifftshift_2D
    import utils.crop_pad


    Nout = outsize;
    Nin = size(im);

    imFT = fftshift_2D(fft2(im));

    imout = crop_pad(imFT, outsize);

    % add noise to avoid correlation between upsampled imaged from
    % interpolation artefacts 
    imout = imout + noise_amplitude * randn(outsize).*mean(min(abs(imFT))); 

    imout = ifft2(ifftshift_2D(imout))*(Nout(1)*Nout(2)/(Nin(1)*Nin(2)));


end

