% REPORT_REFINED_GEOMETRY report results of the geometry refinenement in a readable way 
%
% p = report_refined_geometry(self, param, p)
%
% 
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** param       structure containing parameters for the engines 
% ** p          ptychoshelves p structure 
%
% returns: 
% ** p       updated ptychoshelves p structure 



% Academic License Agreement

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


function p = report_refined_geometry(self, param, p)

    import utils.*
    scale = 1; 
    
   %% GENERATE REPORT ABOUT GEOMETRY REFINEMENT 
   if isempty(p.affine_matrix)
       p.affine_matrix = diag([1,1]); 
   end
   % aux function for printing results 
   mat2str=@(matrix)sprintf(' [%.4g , %.4g ; %.4g , %.4g ] ', reshape(matrix',[],1));

   
   if ~isempty(self.affine_matrix) && param.probe_position_search < param.number_iterations && ~isempty(param.probe_geometry_model)
       for ii = 1:length(self.affine_matrix)
            self.affine_matrix{ii} = rot90(self.affine_matrix{ii},2)';  % rotation is important to match the coordinates with other engines  
       end
       
       for ii = 1:length(self.affine_matrix)
            p.affine_matrix_refined{ii} = p.affine_matrix * self.affine_matrix{ii};
       end        
      if param.Nscans == 2 && param.share_object && param.mirror_objects
        %% use mirrored scans to refine scanning geometry 
        verbose(0, '========================================================= ')
        verbose(0, '==== Geometry parameters for shared 0/180 deg scans ===== ')
        verbose(0, '========================================================= ')
        verbose(0, '')
        % find difference between 0 and 180 , 
        affine_mat_relative = sqrtm(self.affine_matrix{1} * self.affine_matrix{2})*p.affine_matrix;  
        % keep only nondiagonal terms 
        affine_mat_relative = eye(2) + (1-eye(2)).*affine_mat_relative;

        verbose(0, '=============== RELATIVE (0vs180deg) GEOMETRY REFINEMENT ===============')
        verbose(0, '(apply p.affine_matrix manually to your template)')
        verbose(0, 'p.affine_matrix = %s ', mat2str(affine_mat_relative))
        [~, ~, rotation, shear] = math.decompose_affine_matrix(affine_mat_relative); 
        verbose(0, 'This correponds to the following parameters: [rotation=%.3fdeg , shear=%.3fdeg] ', [rotation, shear])

        % find affine matrix that stays contant when moving from 0 to
        % 180 deg, include also the diagonal terms from original affine
        % matrix 
        affine_mat_global = sqrtm(self.affine_matrix{1} * ( [1,-1;-1,1] .* self.affine_matrix{2})); 
        affine_mat_global = affine_mat_global* diag(diag(p.affine_matrix)); 

        verbose(0, '====================================================================================')
        verbose(0, '')
        
        scale = mean(diag(affine_mat_global));
        
      else
        %% use conventional scans to refine scanning geometry 
        median_affine_matrix = median(cat(3,p.affine_matrix_refined{:}),3); 
        verbose(0, '')
        verbose(0, '========= 2D PTYCHO GEOMETRY REFINEMENT, apply manually to your template ===========')
        verbose(0, 'p.affine_matrix = %s' , mat2str(median_affine_matrix))
        verbose(0, '====================================================================================')
        verbose(0, '')
        verbose(0, 'Advanced: ======================== AFFINE CORRECTION OF SCANNER AXIS ====================')
        verbose(0, 'Advanced: (for control system of piezo scanner, important for calibration of cSAXS fast FZP scanner)')
        verbose(0, 'Advanced: correction_matrix = inv(p.affine_matrix) = %s ', mat2str(inv(median_affine_matrix)))
        verbose(0, 'Advanced: ===============================================================================')
        verbose(0, 'Note:  Use scans at 0 and 180 deg with eng.share_object == true && eng.mirror_objects == true to get estimation of the 0vs180deg affine matrix requied for ptychotomography')
        verbose(0, '')
        verbose(0, '')
        verbose(0, '==== Geometry parameters for each scan===== ')
        for ii = 1:length(p.affine_matrix_refined)
           [scale, asymmetry, rotation, shear] = math.decompose_affine_matrix(p.affine_matrix_refined{ii}); 
           verbose(0, 'Scan #%i: [scale=%.4f , asymmetry=%.3f , rotation=%.3fdeg , shear=%.3fdeg, shift = %.1f %.1fpx ] ', [p.scan_number(ii), scale, asymmetry, rotation, shear, self.shift_scans(:,ii)'])
        end  
        scale = mean(diag(median_affine_matrix));
      end

      %% evaluate results if the simulated geometry 
      if isfield(p,'simulation') && check_option(p.simulation,'affine_matrix')
          % report for simulation
        verbose(-2, '')
        verbose(-2, '========== IDEAL AFFINE MATRIX vs RECONSTRUCTED AFFINE MATRIX ====')
        verbose(-2, 'ideal_affine_matrix   = %s ', mat2str(p.simulation.affine_matrix))
        if param.Nscans == 2 && param.share_object && param.mirror_objects
           affine_mat = diag(diag(affine_mat_global)) + affine_mat_relative - eye(2);
        else
           affine_mat = median_affine_matrix;
        end
        verbose(-2, 'refined_affine_matrix = %s ', mat2str(affine_mat))
        verbose(-2, '==================================================================')
        verbose(-2, '')
      end
   end
   
   
   if param.number_iterations > param.detector_rotation_search && ~isempty(param.probe_geometry_model)
         if isfield(p,'simulation') && check_option(p.simulation,'sample_rotation_angles') 
             % report for simulation
            verbose(-2, '')
            verbose(-2, '==== SIMULATION: IDEAL vs RECONSTRUCTED DETECTOR ROTATION CORRECTION  =======')
            verbose(-2, 'ideal camera rotation = %.3f deg     reconstructed camera rotation = %.3f deg', p.simulation.sample_rotation_angles(3), self.detector_rotation(1))
            verbose(-2, '=============================================================================')
            verbose(-2, '')
         else
             % report for real data 
            verbose(0, '')
            verbose(0, '==========  RECONSTRUCTED DETECTOR ROTATION CORRECTION  =====================')
            verbose(0, '(misalignement between detector and the rotation axis, correct by camera rotation)')
            verbose(0, 'Reconstructed camera rotation = %.3f deg', self.detector_rotation(1) + param.sample_rotation_angles(3))
            verbose(0, '=============================================================================')
            verbose(0, '')
         end
   end
   
   
  if param.number_iterations > param.detector_scale_search && ~isempty(param.probe_geometry_model)
         if isfield(p,'simulation') && isfield(p.simulation, 'affine_matrix') && param.detector_scale_search
             % report for simulation
            verbose(-2, '')
            verbose(-2, '============ SIMULATION: IDEAL vs RECONSTRUCTED DETECTOR SCALE  =============')
            if check_option(p.simulation, 'z')
                scale_z = p.z / p.simulation.z; 
            else
                scale_z = 1; 
            end
            verbose(-2, 'ideal scale = %.3f   reconstructed scale = %.3f ', 1/(mean(diag(p.simulation.affine_matrix)) * scale_z), scale/self.detector_scale)
            verbose(-2, '=============================================================================')
            verbose(-2, '')
         else
            % report for real data 
            verbose(0, '')
            verbose(0, '==========  RECONSTRUCTED DETECTOR SCALE CORRECTION  ========================')
            verbose(0, '(relative scaling error of the provided reconstruction pixel p.dx_spec )')
            verbose(0, 'reconstructed scale = %.3f ', scale/self.detector_scale)
            verbose(0, '=============================================================================')
            verbose(0, '')
         end
  end
            
  
 if isinf(self.z_distance) && ...
         ((param.detector_scale_search < param.number_iterations) ...
         || (param.probe_position_search < param.number_iterations && any(ismember(param.probe_geometry_model, 'scale'))))

    verbose(-2, '')
    verbose(-2, '==========  RECONSTRUCTED DETECTOR DISTANCE CORRECTION  =====================')

    if isfield(p,'simulation') && check_option(p.simulation, 'z')
        % report for simulation
        verbose(-2, '==== Compare ideal (simulated) distance and distance refined by ptychography')
        if   isfield(p.simulation, 'affine_matrix') 
             aff_corr_scale = mean(diag(p.simulation.affine_matrix));
        else
             aff_corr_scale = 1; 
        end
        verbose(-2, 'ideal camera distance = %.4f   estimated camera distance = %.4f', p.simulation.z/aff_corr_scale,  p.z / (scale * self.detector_scale))
    else
        % report for measurements 
        verbose(-2, '(needs to be corrected by adjusting p.z parameter in the template)')
        verbose(-2, '==== Scale error corresponds to the following p.z value')
        verbose(-2, 'p.z = %.4f  (error=%.2g%%)',  p.z/(scale*self.detector_scale), 100*(1/(scale*self.detector_scale)-1))
        if param.probe_position_search < param.number_iterations && ~isempty(param.probe_geometry_model)
             verbose(0, '(corrected p.affine_matrix to be used with the new p.z value, add manually to your template)')
             if exist('median_affine_matrix', 'var')
                 affine_mat = median_affine_matrix;
             else
                 affine_mat = affine_mat_relative;
             end
             verbose(0, 'p.affine_matrix = %s ', mat2str(affine_mat / (scale/self.detector_scale) ))
        end
    end
    verbose(-2, '=============================================================================')

 end

end

