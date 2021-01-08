%  ESTIMATE_REQUIRED_GPU_MEMORY Estimate GPU memory required to run reconstruction with provided parameters 
%
% [required_mem , data_mem, object_mem, required_fft_mem] = ...
%                   estimate_required_GPU_memory(self, par)
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** par       structure containing parameters for the engines 
% returns 
% ++  required_mem           total required mem
% ++  data_mem               mem to store data 
% ++  object_mem             mem to store object 
% ++  required_fft_mem       mem to run FFT 


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


function [required_mem , data_mem, object_mem, required_fft_mem] =  estimate_req_memory(self, par, grouping)
        % estimate how much GPU memoty will be needed to run the code 
        import utils.*
        import engines.GPU_MS.shared.*
        if nargin < 3
            grouping = par.grouping; 
        end
        if par.gpu_id < 1
            % ie no GPU is used 
            required_mem = nan; data_mem = nan; object_mem = nan; required_fft_mem = nan; 
            return
        end
        
        
        required_mem = 0; 
        data_mem = 0;
        object_mem = 0; 
        verbose(2, 'Checking available GPU memory ID:%i', par.gpu_id)
        if par.keep_on_gpu
            data_class = class(self.diffraction); 
            switch data_class
                case 'single'
                    byte_size = 4;
                case 'uint16'
                    byte_size = 2;
                case {'uint8','int8'}
                    byte_size = 1;
            end
            % keep data on GPU
            data_mem = data_mem + prod(self.Np_p)*self.Npos*byte_size / 2^(2*par.upsampling_data_factor); % self.diffraction can be either cell or array, so the size is calculated from par
            data_mem = data_mem + numel(self.mask);  % bool (uint8 in matlab)
            data_mem = data_mem + numel(self.noise)*4; % single 
            % keep views on GPU
            if any(strcmpi(par.method, {'DM'}))
                required_mem = required_mem + 8*prod(self.Np_p)*self.Npos*par.probe_modes;
            end
        end
        % very empirical guess , assuming FFT memory requirement ~6*8*numel(x)
        required_fft_mem = (2)*6*8*prod(self.Np_p)*grouping;  % empirically tested  
        required_mem = required_mem  + required_fft_mem; 

        % basic arrays: obj_proj, chi
        % add for multiple layers by ZC
        if isfield(par,'Nlayers')
            Nlayers=par.Nlayers;
        else
            Nlayers=1;
        end
        
        required_mem = required_mem + 2*4*2*prod(self.Np_p)*grouping*par.Nmodes; % added Nmodes by ZC
                
        % multilayers of probe added by ZC
        required_mem = required_mem + 2*4*prod(self.Np_p)*Nlayers*grouping*par.Nmodes + 2*4*prod(self.Np_p)*(grouping+par.Nmodes-1);
        
        if any(~isinf([par.probe_position_search, par.probe_fourier_shift_search]))
           % memory needed to keep a probe for each position in the grouping 
           required_mem = required_mem + 2*4*prod(self.Np_p)*grouping;
        end
        if is_method(par, 'PIE') && par.variable_probe
            % size of probe stored for each scan position 
            required_mem = required_mem + 4*2*prod(self.Np_p)*self.Npos;
        end

        if is_method(par, 'ML') && par.momentum
            % account for data stored for momentum estimate 
            momentum_mem = 2 + 1 ; % 2 previous steps + velocity map
            object_mem = object_mem + momentum_mem*4*2*prod(self.Np_o)*numel(self.object);
        end
        if is_method(par, 'ML') && par.accelerated_gradients_start < par.number_iterations
            % account for data stored for accelerated gradient
            object_mem = object_mem + 2*4*2*prod(self.Np_o)*numel(self.object);
        end
        if ~par.share_probe
            % size of probe stored for each scan position 
            required_mem = required_mem + 4*2*prod(self.Np_p)*grouping;
        end
        if par.apply_subpix_shift || par.variable_probe
            % memory needed for subpixel shifted probe 
            required_mem = required_mem + 4*2*prod(self.Np_p)*grouping*par.Nmodes; % add par.Nmodes by ZC
        end
        % size of the object, object update, local object illumination,
        % total object illumination 
        object_mem = object_mem + 4*prod(self.Np_o) * ( (2+2+1)*numel(self.object)  ); 
        
        % add object and data 
        required_mem = required_mem + data_mem;
        required_mem = required_mem+object_mem;
end
