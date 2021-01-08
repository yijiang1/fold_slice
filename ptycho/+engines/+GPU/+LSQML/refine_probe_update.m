% REFINE_PROBE_UPDATE calculate improved update direction
% apply "overlap" constraint to get better estimate  of the update direction
% also updates the variable probe estimate 
%
% [ self,m_probe_update, probe_update, cache]  = ...
%                 refine_probe_update(self, obj_proj, probe_update, chi,layer_ids,probe_id,p_ind,g_ind, par, cache)
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** obj_proj  [Nx, Ny, N] array, estimate of the object update for each scan position, ie conj(P)*chi
% ** probe_update      cell of object sized arrays containg previous optimal updates 
% ** chi      indices corresponding to the current group that is solved in parallel 
% ** layer_ids  id of the solved layer for multilayer ptycho 
% ** probe_id   id of incoherent probe mode 
% ** p_ind      indices containg corresponding probe id for each processed position
% ** g_ind      indices corresponding to the current group that is solved in parallel 
% ** par       structure containing parameters for the engines 
% ** cache     structure with precalculated values to avoid unnecessary overhead
%
% returns:
% ++ object_upd_sum      cell of object sized arrays containg updated optimal update
% ++ object_update_proj  [Nx, Ny, N] array, estimate of the refiend object update for each scan position,
% ++ cache                structure with precalculated values to avoid unnecessary overhead
%
%
% see also: engines.GPU.LSQML 


function [ self,m_probe_update, probe_update, cache]  = ...
                refine_probe_update(self, obj_proj, probe_update, chi,layer_ids,probe_id,p_ind,g_ind, par, cache)
    % get probe update direction 
    import engines.GPU.shared.*
    import engines.GPU.LSQML.*
    import engines.GPU.GPU_wrapper.*
    import math.*
    import utils.*

    if layer_ids > 1 % in case of multilayer object 
        m_probe_update = []; 
        return
    end
    
    if (probe_id == 1 && par.apply_subpix_shift && isinf(self.z_distance(end))) || is_used(par,'fly_scan')
        probe_update  = apply_subpx_shift(probe_update , -self.modes{min(end,probe_id)}.sub_px_shift(g_ind,:) );
    end
    if probe_id == 1
        probe_update = apply_subpx_shift_fft(probe_update, -self.modes{min(end,probe_id)}.probe_fourier_shift(g_ind,:));
    end

  
    if par.share_probe || length(unique(p_ind)) == 1 
       % BETTER WAY: assume that sum(|obj_proj|^2,3) is close to 1 
       % and additionally use weighting based on confidence given by illum_sum_0
        % => apriory weighting giving less importance to the less
        % illuminated regions
        %weight_proj = cache.illum_sum_0{1} ./ (cache.illum_sum_0{1}+0.01*cache.MAX_ILLUM(1));
        %weight_proj = get_views({weight_proj},[],1,1, g_ind, cache);
        %m_probe_update = mean( weight_proj.* probe_update,3);
        % or originally was used simple average , good for object >> probe 
        m_probe_update = mean(probe_update,3);  % calculate single update for all current positions
    else  % unshared probe and multiple scans in one group (ie shared object)
        probe_ids = unique(p_ind); 
        m_probe_update = Gzeros([self.Np_p,length(probe_ids)], true);
        for probe_id =probe_ids(:)'
            m_probe_update(:,:,probe_id) = mean(probe_update(:,:,p_ind == probe_id),3); % calculate one update for each scan
        end
    end

    if (par.variable_probe || par.variable_intensity) && probe_id == 1  
        % ORTHOGONAL PROBE RELAXATION (OPRP) EXTENSION - allow
        % variable probe wavefront 
        %  Odstrcil, M., et al. "Ptychographic coherent diffractive imaging with orthogonal probe relaxation." Optics express 24.8 (2016): 8360-8369.
        % iterate over all sub probes 
        [self.probe{probe_id}, self.probe_evolution] = ...
            update_variable_probe(self, self.probe{probe_id}, self.probe_evolution, m_probe_update, probe_update, obj_proj, chi,cache.illum_sum_0{probe_id}, p_ind, g_ind, cache, par);
    end

%     % apply probe constraints 
%     if probe_id == 1 && (check_option(self,'probe_support') || check_option(self,'probe_support_fft'))
%         m_probe_update = apply_probe_contraints(m_probe_update, self.modes{probe_id});
%     end

end