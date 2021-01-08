% SAVE_TO_P save parameters and recosntrutions from  param and self structures to the p-structure
%
% p_out = save_to_p(self, param, p, fourier_error)
%
% 
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** param       structure containing parameters for the engines 
% ** p          ptychoshelves p structure 
% ** fourier_error  array [Npos,1] containing evolution of reconstruction error 
%
% returns: 
% ** p_out       updated ptychoshelves p structure 



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


function p_out = save_to_p(self, param, p, fourier_error)
   
    import utils.*
    Np_p = self.Np_p; 

    p_out = p; 
    

    %% return calculated values back to the main code (p-structure)
    p_out.probe_modes = length(self.probe);
    p_out.numprobs = size(self.probe{1},3);
    for ll = 1:p_out.probe_modes
        if size(self.probe{ll},3) == self.Npos
            % classical OPRP method 
            for ii = 1:p.numscans
                probes(:,:,ii,ll) = mean(self.probe{ll}(:,:,self.reconstruct_ind{ii}),3);
            end
             % save the variable modes 
             [X,V] = extract_PCA(self.probe{ll}, param.variable_probe_modes);
             p_out.probe_PCA.eigen_vec = X * (prod(sqrt(Np_p))*2*p.renorm);% normalization for consistency with the CPU code
             p_out.probe_PCA.evolution = V;
        else
            % store constant part 
            probes(:,:,1:size(self.probe{ll},3),ll) = self.probe{ll}(:,:,:,1);
            if ndims(self.probe{ll}) == 4
                p_ind = zeros(self.Npos,1);
                for kk = 1:length(self.reconstruct_ind)
                   p_ind = p_ind + kk*ismember(1:self.Npos,  self.reconstruct_ind{kk})';
                end
                p_out.probe_variable.eigen_vec = self.probe{ll}(:,:,:,2) * (prod(sqrt(Np_p))*2*p.renorm);% normalization for consistency with the CPU code
                p_out.probe_variable.evolution = self.probe_evolution;
            end
        end
    end
            
    
    % normalization for consistency with the CPU code
    probes = probes * (prod(sqrt(Np_p))*2*p.renorm); 
    p_out.probes = probes; 
   
    
    position_offset = 1+floor((p.object_size-self.Np_p)/2);

    % for consistency with the CPU code revert the object to the original size 
    p_out.numobjs = size(self.object,1);
    p_out.object = cell(1,p_out.numobjs);
    for i = 1:p_out.numobjs
        obj_size = p.object_size(min(end,p.share_object_ID(i)),:); 
        p_out.object{i} = single([]); 
        for layer = 1:param.Nlayers
            p_out.object{i}(:,:,1,layer) = imshift_fast(self.object{i,layer}, -1,-1, obj_size, 'nearest', mean(self.object{i,layer}(:)));
        end
    end

    for i = 1:p.numscans
        id = p.share_object_ID(i); 
        obj_size = p.object_size(min(end,id),:);
        p_out.illum_sum{id} = imshift_fast(self.illum_sum{id}, -1,-1, obj_size, 'nearest');
    end
    
    if param.probe_position_search < param.number_iterations
        % store the refined positions 
        p_out.positions = self.probe_positions(:,[2,1]);
        % return to the original coordinates 
        p_out.positions_0 = self.probe_positions_0(:,[2,1]);
        for i = 1:length(self.reconstruct_ind)
            ind = self.reconstruct_ind{i};
            p_out.positions(ind,:) = p_out.positions(ind,:) + position_offset(p.share_object_ID(i),:);
            p_out.positions_0(ind,:) = p_out.positions_0(ind,:) + position_offset(p.share_object_ID(i),:);
        end
    else
        for i = 1:length(self.reconstruct_ind)
            ind = self.reconstruct_ind{i};
            if ~isempty(self.probe_positions)
                p_out.positions(ind,:) = self.probe_positions(ind,[2,1]) + position_offset(p.share_object_ID(i),:);
            else
                p_out.positions(ind,:) = self.probe_positions_0(ind,[2,1]) + position_offset(p.share_object_ID(i),:);
            end
        end
    end
    

    if param.probe_position_search < param.number_iterations || param.detector_rotation_search < param.number_iterations  || param.detector_scale_search < param.number_iterations
        p_out = engines.GPU.analysis.report_refined_geometry(self, param, p_out); 
    end
 

    % save additional reconstructed parameters  
    for item = {'background', 'intensity_corr', 'probe_fourier_shift' }
        try
            p_out.(item{1}) = self.(item{1});
        end
    end
    
    % save error metrics 
    ind_ok = any(~isnan(fourier_error),2); % plot only the reported values 
    p_out.error_metric.value = nanmean(fourier_error(ind_ok,:),2);
    p_out.error_metric.iteration = find(ind_ok);
    if strcmp(param.likelihood,'poisson' )
        p_out.error_metric.err_metric = 'poisson';
    else
        p_out.error_metric.err_metric = 'L1';
    end
    p_out.error_metric.method = ['GPU-',param.method, ' metric:' , param.likelihood ];

    
end
