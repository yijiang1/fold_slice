% GRADIENT_NF_PROPAGATION_SOLVER update estimate of the nearfield propagation distance based on the 
% current update of the wavefront 
%
% self = gradient_NF_propagation_solver(self,psi,chi, cache, ind)
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** psi       exitwave (psi = P*O)
% ** chi       [Nx,Ny,N] array, difference between original and updated exit-wave  
% ** cache     precalculated values
% ** ind      indices containg corresponding probe id for each processed position
%
% returns:
% ++ self    updated self structure with optimized propagation distance 
% ++ cache    precalculated values with updated velocity for momentum method 


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
%   



function [self, cache] = gradient_NF_propagation_solver(self,psi,chi, cache, ind, layer_id)
    % update estimate of the nearfield propagation distance based on the 
    % current update of the wavefront 

    import engines.GPU_MS.GPU_wrapper.*
    import engines.GPU_MS.shared.*
    import math.*
    import utils.*

    dH = cache.ASM_difference;
    
    % propagate to farfield 
    Psi_0 = fft2_safe(psi{1});
    Psi_1 = fft2_safe(psi{1}+chi{1});
    % get phase diffence 
    Psi_diff = sum(Psi_0 .* conj(Psi_1),3);
    aPsi = abs(Psi_diff);
    Psi_diff = Psi_diff ./ (aPsi+mean2(aPsi).*1e-6);
    
    % calculate distance that best explains the differences 
    dz = -Ggather(sum2(real(conj(dH) .* Psi_diff)) ./ sum2(abs(dH).^2));
    
     
%     %% USE MOMENTUM ACCELERATION TO MAKE THE CONVERGENCE FASTER
%     try
%     momentum_memory = 10; % compare 10 iterations 
%     if length(self.modes{1}.distances) > momentum_memory
%         if ~isfield(cache, 'velocity_NF_propagation')
%             cache.velocity_NF_propagation = 0; 
%         end
%         
%         updates = diff(self.modes{1}.distances(end-momentum_memory:end)); 
%         corr_level = corr(updates(1:end-1)', updates(2:end)'); 
%         
%         corr_level
%         
%         
%         if all(corr_level > 0 )
%             %%%%%%%%%%%% very empritical model %%%
%             gain = 1;                           % smaller -> lower relative speed (less momentum)
%             friction =  0.1*max(0, 0.5-corr_level);   % smaller -> longer memory, more momentum 
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         else
%            gain = 0; friction = 0.5; 
%         end
%         
%         cache.velocity_NF_propagation = cache.velocity_NF_propagation*(1-friction) + dz;
%         %% apply the velocity to the refined positions , if the postition updated are sufficiently small 
% 
%         cache.velocity_NF_propagation
% 
%         dz = dz + gain*cache.velocity_NF_propagation; 
%         
% 
%     end
%     catch
%         keyboard
%     end
    
    
    % update propagation values 
    if any(ind==1)
        % add new field only when new iteration is started 
        self.modes{1}.distances = [self.modes{1}.distances,self.modes{1}.distances(end) + dz];
        verbose(1, 'Propagation distance: %3.5gum ', self.modes{1}.distances(end)*1e6)
    else
        self.modes{1}.distances = self.modes{1}.distances + dz;
    end
    
    % use the same distance for all positions 
    for i = 1:length(self.modes)  
        self.modes{i}.distances = self.modes{1}.distances;
    end

end