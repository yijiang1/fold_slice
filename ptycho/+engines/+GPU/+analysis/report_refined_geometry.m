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
            %switch diagonal elements
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

