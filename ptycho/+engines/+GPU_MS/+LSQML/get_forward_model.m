% GET_FORWARD_MODEL from the provided object and probe calculate the exit wave 
%
% [self, probe, obj_proj, psi] = get_forward_model(self, obj_proj, par, cache, g_ind, p_ind, scan_ids, layer_ids)
% 
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** obj_proj  [Nx,Ny,N] array, just a preallocated array on GPU, can be empty 
% ** par       structure containing parameters for the engines 
% ** cache     structure with precalculated values to avoid unnecessary overhead
% ** g_ind      indices corresponding to the current group that is solved in parallel 
% ** p_ind      indices containg corresponding probe id for each processed position
% ** scan_ids   determines to which scan correponds each of the position 
% ** layer_ids  id of the solved layer for multilayer ptycho 
%
%
% returns:
% ++ self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ++ probe     either [Nx,Nx,1] or [Nx,Nx,N] aarray of shared probe or variable probe that differs for each position 
% ++ obj_proj  [Nx,Ny,N] array, views of the object for each scan position 
% ++ psi       [Nx,Ny,N] array, complex valued exit-wave (psi = P*O)
%
% see also: engines.GPU_MS.LSQML 

function [self, probe, obj_proj, psi] = get_forward_model(self, obj_proj, par, cache, g_ind, p_ind, scan_ids, layer_ids)
    import engines.GPU_MS.shared.*
    import engines.GPU_MS.GPU_wrapper.*
    import engines.GPU_MS.LSQML.*
    import math.*
    import utils.*
    import plotting.*

    if isempty(obj_proj{1})
        for ll = 1:par.object_modes
            obj_proj{ll} = Gzeros([self.Np_p, 0], true);
        end
    end
    
    % allocate memory for probes. Added by ZC
    % probe = self.probe;
    if par.Nlayers > 1
        probe=cell(par.probe_modes,par.Nlayers+1);
    end
    
    % get illumination probe 
    for ll = 1:par.probe_modes
        if (ll == 1 && (par.variable_probe || par.variable_intensity))
            % add variable probe (OPRP) part into the constant illumination 
            probe{ll,1} =  get_variable_probe(self.probe{ll}, self.probe_evolution(g_ind,:),p_ind{ll});
        else
            % store the normal (constant) probe(s)
            probe{ll,1} = self.probe{min(ll,end)}(:,:,min(end,p_ind{ll}),1);
        end

        if (ll == 1 && par.apply_subpix_shift && isinf(self.z_distance(end)))  || is_used(par,'fly_scan')
            % only in farfield mode 
            probe{ll,1} = apply_subpx_shift(probe{ll,1}, self.modes{min(end,ll)}.sub_px_shift(g_ind,:) );
        end
        if (ll == 1)
            probe{ll,1} = apply_subpx_shift_fft(probe{ll,1}, self.modes{1}.probe_fourier_shift(g_ind,:)); 
        end
    end
    
   % allocate memory for psi. Added by ZC
   psi=cell(max(par.object_modes, par.probe_modes),1);

   % get projection of the object and probe 
   for layer = 1:par.Nlayers
       for ll = 1:max(par.object_modes, par.probe_modes)
           llo = min(ll, par.object_modes); 
           llp = min(ll, par.probe_modes); 
            % get objects projections 
            obj_proj{llo} = get_views(self.object, obj_proj{llo},layer_ids(layer),llo, g_ind, cache, scan_ids,[]);
            if (ll == 1 && par.apply_subpix_shift && ~isinf(self.z_distance(end)))
                % only in nearfield mode , apply shift in the opposite direction 
                obj_proj{ll} = apply_subpx_shift(obj_proj{ll} .* cache.apodwin, -self.modes{min(end,ll)}.sub_px_shift(g_ind,:) ) ./ cache.apodwin;
            end
            
            % get exitwave after each layer
            psi{ll} = probe{llp,layer} .* obj_proj{llo};
            
            %modified by YJ: no need to do another propagation after the
            %last object layer
            if layer < par.Nlayers 
                % fourier propagation
                [psi{ll}] = fwd_fourier_proj(psi{ll} , self.modes{layer}, g_ind);  
                if par.Nlayers > 1
                     probe{llp,layer+1} = psi{llp};
                end
            end
       end
   end
   % At this point, psi is the exit wave after the last object layer
   % Now propagate the wave function to far-field (detector) plane
   % This is same as using fwd_fourier_proj w. distance = inf and no camera angle refinement 
   for ll= 1:max(par.object_modes, par.probe_modes)
       psi{ll} = fft2_safe(psi{ll});  % fully farfield 
   end
end
