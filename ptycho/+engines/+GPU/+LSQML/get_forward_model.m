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
% see also: engines.GPU.LSQML 


function [self, probe, obj_proj, psi] = get_forward_model(self, obj_proj, par, cache, g_ind, p_ind, scan_ids, layer_ids)
    import engines.GPU.shared.*
    import engines.GPU.GPU_wrapper.*
    import engines.GPU.LSQML.*
    import math.*
    import utils.*
    import plotting.*

    if isempty(obj_proj{1})
        for ll = 1:par.object_modes
            obj_proj{ll} = Gzeros([self.Np_p, 0], true);
        end
    end
    probe = self.probe;
    %length(probe)=par.probe_modes
    %size(probe{1})=[Np,Np,1,par.variable_probe_modes+1];
    
    % get illumination probe 
    for ll = 1:par.probe_modes
        %p_ind{ll} is always 1 for single dataset
        if (ll == 1 && (par.variable_probe || par.variable_intensity))
            % add variable probe (OPRP) part into the constant illumination
            % OPRP only applies to the FIRST probe mode in mixed-states
            probe{ll,1} =  get_variable_probe(self.probe{ll}, self.probe_evolution(g_ind,:),p_ind{ll});
            % size(probe{ll,1}) = [Np_p(1), Np_p(2), # of probes in this group];
        else
            % store the normal (constant) probe(s)
            probe{ll,1} = self.probe{min(ll,end)}(:,:,min(end,p_ind{ll}),1);
            % size(probe{ll,1}) = [Np_p(1), Np_p(2), 1]; No OPR for higher probe modes
        end
    
        if (ll == 1 && par.apply_subpix_shift && isinf(self.z_distance(end))) || is_used(par,'fly_scan')
            % only in farfield mode 
            probe{ll} = apply_subpx_shift(probe{ll}, self.modes{min(end,ll)}.sub_px_shift(g_ind,:) );
        end
        if (ll == 1)
            probe{ll} = apply_subpx_shift_fft(probe{ll}, self.modes{1}.probe_fourier_shift(g_ind,:)); 
        end 
    end

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
            % fourier propagation  
            [psi{ll}] = fwd_fourier_proj(psi{ll} , self.modes{layer}, g_ind);  
            if par.Nlayers > 1
                 probe{llp,layer+1} = psi{llp};
            end
       end
   end
end
