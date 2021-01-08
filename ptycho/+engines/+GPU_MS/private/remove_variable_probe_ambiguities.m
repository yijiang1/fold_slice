% REMOVE_VARIABLE_PROBE_AMBIGUITIES remove ambiguities in when the variable probe methods is used 
% -> normalize the OPR modes 
% -> make sue that the modes are orthogonal 
% 
% self = remove_variable_probe_ambiguities(self,par)
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** par       structure containing parameters for the engines 
%
% returns:
% ++ self        self-like structure with final reconstruction
%
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


function self = remove_variable_probe_ambiguities(self,par)
    % apply normalization and orthogonalize the variable probe modes 
    import math.*
    import engines.GPU_MS.GPU_wrapper.*
    Nmodes = par.variable_probe_modes; 
    if par.share_probe
        inds = {[self.reconstruct_ind{:}]};
    else
        inds = self.reconstruct_ind;
    end
    probe_modes = length(inds);   % number of independend probes 
    probe = self.probe{1}; 
    
    for kk = 1:probe_modes
        ind = inds{kk};
        if par.variable_probe_smooth > 0
            % relaxed regularization , restrict the probe evolution to slow
            % changes with polynomial order given by par.variable_probe_smooth
            for ll = 2:Nmodes
                self.probe_evolution(ind,ll) = self.probe_evolution(ind,ll) * 0.5 + 0.5*polyval(polyfit(ind,self.probe_evolution(ind,ll)',round(par.variable_probe_smooth)), ind)';
            end
        end
        % remove degrees of freedom (orthogonalize the first mode)
%         self.probe_evolution(ind,:) = self.probe_evolution(ind,:) - mean(self.probe_evolution(ind,:),1); 
%         self.probe_evolution(ind,1) = self.probe_evolution(ind,1)*0.99 + 1; 
    end

    %% apply normalization on the coherence modes 
    vprobe_norm = Ggather(norm2(probe(:,:,:,2:end)));
    probe(:,:,:,2:end) = probe(:,:,:,2:end) ./ vprobe_norm;
    for kk = 1:probe_modes
        ind = inds{kk};
        self.probe_evolution(ind,2:end) = self.probe_evolution(ind,2:end) .* reshape(vprobe_norm(1,1,kk,:),1,[]);
    end
    
    %% orthogonalize the variable modes
    for i = 1:(1+Nmodes)
       for j = 1:i-1
           mx_j = probe(:,:,:,j);
           mx_i = probe(:,:,:,i);
           proj = sum2(mx_i .* conj(mx_j)) ...
                ./ sum2(abs(mx_j).^2); 
           probe(:,:,:,i)  =  probe(:,:,:,i) -  proj .* mx_j;
           
           for kk = 1:probe_modes
               ind = inds{kk};
               p_j = self.probe_evolution(ind,j);
               p_i = self.probe_evolution(ind,i);
               proj = sum2(p_i .* conj(p_j)) ...
                    ./ sum2(abs(p_j).^2); 
               self.probe_evolution(ind,i)  =  self.probe_evolution(ind,i) -  proj .* p_j;
           end
       end
    end
    

   %% sort the modes by their energy 
   Energy = zeros(probe_modes, Nmodes); 
   for kk = 1:probe_modes
       for ii = 1:Nmodes
            ind = inds{kk};
            Energy(kk,ii) = norm(self.probe_evolution(ind,ii+1));
       end
   end
   [~,ind] = sort(Energy,2, 'descend'); 
   for kk = 1:probe_modes
       swap = [1,1+ind(kk,:)]; 
       probe(:,:,kk,:) = probe(:,:,kk,swap); 
       self.probe_evolution(inds{kk},:) = self.probe_evolution(inds{kk},swap); 
   end

   %% remove outliers 
   aevol = abs(self.probe_evolution); 
   self.probe_evolution = min(aevol, 1.5*quantile(aevol,0.95)) .* sign(self.probe_evolution); 

    % store the results 
    self.probe{1} = probe; 
    
end


