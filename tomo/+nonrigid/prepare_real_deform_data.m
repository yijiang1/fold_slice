% prepare_real_deform_data - load already prealigned and fixed projections
%
%[dphase, shift_3D_all_0, angles, par] = ...
%                 prepare_real_deform_data(path_to_projections, projection_filename, sample_name, par_0)
%
% Inputs:
%    **path_to_projections     path where are stored preloaded data  
%    **projection_filename     name of the file where are the preloaded data 
%    **sample_name             name of the sample, it used for cache file reparation 
%    **par_0                   initial paramters structure that will be merged with the loaded paramters 
% Outputs: 
%    ++dphase               phase difference for the complex project 
%    ++shift_3D_all_0       original deformation = {}
%    ++angles               angles for each of the projection 
%    ++par                  merged parameter structure 
%    ++object               complex valued projections 

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2018 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%*-----------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing language this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS matlab package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.

function [dphase, shift_3D_all_0, angles, par, object] = ...
                 prepare_real_deform_data(path_to_projections, projection_filename, sample_name, par_0)

     import utils.*
     import math.*
     import plotting.*
             
    shift_3D_all_0 = {};
    if ~exist('cache', 'dir')
        mkdir('cache')
    end
    
    cached_file = fullfile(path_to_projections, ['cache_nonrigid_tomo', sample_name,'.mat']); 
    % try to load cached data 

        verbose(0, 'Loading prepared data from %s', fullfile(path_to_projections, projection_filename))
        %% load complex valued projections  
        d = load(fullfile(path_to_projections, projection_filename));
        verbose(0, 'Loading done')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% create data and geometry
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%

         if isfield(d, 'dphase')
            dphase = d.dphase; 
         else
            object = complex(d.stack_object_r,d.stack_object_i); 

            object_ROI = d.object_ROI;

            
            object = object(object_ROI{:},:);
            

            object = smooth_edges(object); 

            dphase = tomo.block_fun(@math.get_phase_gradient_1D,object, 2);

         end
         if isfield(d, 'angles')
            angles = d.angles; 
         else
            angles = d.theta; 
         end
         
         
        dphase = smooth_edges(dphase); 

        
        

        try
            par = d.par; 
        catch
            par = struct(); 
        end
        
        

    % rewrite loaded params by some defaults 
   for field = fields(par_0)'
        par.(field{1}) = par_0.(field{1}); 
   end
   
   par.output_folder = path_to_projections; 
   
   if debug()
        % plot angular blocks 
        figure
        plot(angles, '.-')
        for ii = 1:Nblocks
           plotting.vline(ii*Bsize) 
        end
        xlabel('Projection #')
        ylabel('Angle [deg]')
        title('Angular block splitting')
        drawnow 
   end

   
end