% PIE - generalized version of the ptychographic iterative engine 
%
%[self, cache, fourier_error] =  PIE(self,par,cache,fourier_error,iter)
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** par       structure containing parameters for the engines 
% ** cache     structure with precalculated values to avoid unnecessary overhead
% ** fourier_error  array [Npos,1] containing evolution of reconstruction error 
% ** iter      current iteration number 
% returns:
% ++ self        self-like structure with final reconstruction
% ++ cache     structure with precalculated values to avoid unnecessary overhead
% ++ fourier_error  array [Npos,1] containing evolution of reconstruction error 
%
%
%   Publications most relevant to the Difference-Map implementation
%       Odstrcil, M., Baksh, P., Boden, S. A., Card, R., Chad, J. E., Frey, J. G., & Brocklesby, W. S
%       "Ptychographic coherent diffractive imaging with orthogonal probe relaxation."
%       Optics express 24.8 (2016): 8360-8369.


function [self, cache, fourier_error] =  PIE(self,par,cache,fourier_error,iter)
        import engines.GPU.GPU_wrapper.*
        import math.*
        import utils.*
        import engines.GPU.shared.*
        import engines.GPU.PIE.*
        
        %% MULTILAYER EXTENSION fix or remove in future
        par.probe_modes = max(par.Nlayers, par.probe_modes); 
        par.object_modes = max(par.Nlayers, par.object_modes); 
        
        %disp(size(self.probe{1}))
        if par.variable_probe
            % expand the probe to the full size
            probe_0 = self.probe{1}; 
            self.probe{1} = reshape(self.probe{1},prod(self.Np_p),[]);
            self.probe{1} = reshape(self.probe{1} * self.probe_evolution', self.Np_p(1), self.Np_p(2), []);
        end
        %disp(size(self.probe_evolution))

        par.multilayer_object = par.Nlayers > 1; 
        par.multilayer_probe = false;  % not supported anymore 
        probe_amp_corr = [0,0]; 

        psi = cell(par.Nmodes, 1);
        Psi = cell(par.Nmodes, 1);
        probe_max = cell(par.probe_modes,1); 
        object_max = cell(par.object_modes,1);
        aprobe2 = abs(mean(self.probe{1},3)).^2; % update only once per iteration 
                
        if (is_method(par, 'ePIE') && ...
                ... % empirical estimation when the hybrid PIE method should be used 
                par.grouping > self.Npos/sqrt( pi^2 * mean(Ggather(cache.MAX_ILLUM)) / max(aprobe2(:)))) %  ||  ...% use hybrid ePIE in case of large grouping
%                 (~isempty(self.modes{end}.ASM_factor) && par.grouping > 1) || ...
%                 (is_method(par, 'ePIE')&&  par.multilayer_object)
            par.method = 'hPIE';  % hybrid ePIE
            if iter == 1;verbose(1,'Switching to hybrid PIE method '); end 
        end
        if is_method(par, {'hPIE'})
            for ll = 1:(par.object_modes*par.Nscans)
                for layer = 1:par.Nlayers
                    obj_illum_sum{ll,layer} = Gzeros(self.Np_o);
                    object_upd_sum{ll,layer} = Gzeros(self.Np_o, true);
                end
            end
            for ll = 1:par.probe_modes
                probe_upd_sum{ll}= (1+1i)*1e-8; 
                probe_illum_sum{ll} = 1e-8; 
            end
        end

        
        for ll = 1:par.Nlayers
            obj_proj{ll} = Gzeros([self.Np_p, par.grouping], true);
        end    
        if par.delta_p 
            grad = @get_grad_lsq;  % dumped least squares gradient (preconditioner)
        else
            grad = @get_grad_flat; % common PIE gradient 
        end

        %%%%%%%%%%%%%%%%%% ePIE algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        if par.grouping > 1
            for ll = 1:par.probe_modes
                if ~par.multilayer_object || ll == 1  % update only first layer probe 
                    probe_max{ll} = Ggather(max2(abs(self.probe{ll}))).^2;
                end
            end
            for ll = 1:par.object_modes
                if ~par.multilayer_probe || ll == 1  % update only first layer object 
                    object_max{ll} = Ggather(max2(abs(self.object{ll}(cache.object_ROI{:})))).^2;
                end
            end
        end
        
        if any(isnan(object_max{1})) ||  any(isnan(probe_max{1}))            
            error('Object or probe is nan, try smaller grouping')
        end
    

    
        % use already precalculated indices 
        rand_ind = randi(length( cache.preloaded_indices_sparse));
        indices = cache.preloaded_indices_sparse{rand_ind}.indices;
        scan_ids = cache.preloaded_indices_sparse{rand_ind}.scan_ids;

        for  ind_ii = randperm(length(indices))
            g_ind = indices{ind_ii}; 
            for ll = 1:max([par.probe_modes, par.object_modes,par.Nlayers]) 
                % generate indices of the used probes 
                if par.variable_probe && ll == 1
                    p_ind{ll} = g_ind;
%                 elseif par.multilayer_object && ll > 1
%                     p_ind{ll} = 1:length(ii);
                else  % single probe only
                     if par.share_probe %|| ll > 1  % share incoherent modes 
                        p_ind{ll} = 1;
                     else
                         if all(scan_ids{ind_ii} == scan_ids{ind_ii}(1))
                            p_ind{ll} =  scan_ids{ind_ii}(1);
                         else
                            p_ind{ll} =  scan_ids{ind_ii};
                         end
                     end
                end
            end

            %% load data to GPU (if not loaded yet)
            modF = get_modulus(self, cache, g_ind);
            mask = get_mask(self, cache.mask_indices, g_ind);
            noise = get_noise(self, par, g_ind); 
            

            % get objects projections 
            
            
            
            for layer = 1:par.Nlayers
                ll = 1; 
                obj_proj{layer} = get_views(self.object, obj_proj{layer},layer,ll,g_ind, cache, scan_ids{ind_ii},[]);
            end
            % get illumination probe 
            for ll = 1:par.probe_modes
                if ~par.multilayer_object || ismember(ll, [1,par.Nlayers:par.probe_modes])
                    probe{ll} = self.probe{min(ll,end)}(:,:,min(end,p_ind{ll}));
                end
            end
                       
            for ll = 1:max([par.probe_modes, par.object_modes, par.Nlayers])  % Nlayers 
                %% fourier propagation  
               
                if (ll == 1 && (par.multilayer_object || par.multilayer_probe) )
                    probe{ll} = self.probe{min(ll,end)}(:,:,min(end,p_ind{ll}));
                end
                if ll > 1 && par.multilayer_object ||  par.grouping == 1 % raw ePIE 
                    probe_max{ll} = max(Ggather(max2(abs(probe{ll})))).^2; 
                end
                %disabled by YJ. seems like a bug
                %{ 
                if ll > 1 && par.multilayer_probe ||  par.grouping == 1 % raw ePIE 
                    object_max{ll} = max(Ggather(max2(abs(obj_proj{ll})))).^2;
                end
                %}
                if (ll == 1 && par.apply_subpix_shift) 
                    probe{ll}  = apply_subpx_shift(probe{ll}, self.modes{min(end,ll)}.sub_px_shift(g_ind,:) );
                end
                
                probe{ll} = apply_subpx_shift_fft(probe{ll}, self.modes{1}.probe_fourier_shift(g_ind,:)); 


                % get projection of the object and probe 
                psi{ll} = bsxfun(@times, obj_proj{min(ll,end)}, probe{min(ll,end)});
                Psi{ll} = fwd_fourier_proj(psi{ll} , self.modes{min(end, ll)}); 
                
                if ll < par.Nlayers && par.multilayer_object
                    probe{ll+1} =  Psi{ll};
                end
                
                if ll < par.Nlayers && par.multilayer_probe
                    obj_proj{ll+1} =  Psi{ll};
                end
            end
                        
            % get intensity (modulus) on detector including different corrections
            aPsi = get_reciprocal_model(self, Psi, modF, mask,iter, g_ind, par, cache);
                
                              
            %%%%%%%%%%%%%%%%%%%%%%% LINEAR MODEL CORRECTIONS END %%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%            
            if iter > 0 && par.get_error && (mod(iter,min(20, 2^(floor(2+iter/50)))) == 0 || (iter < 20)) || iter == par.number_iterations  % calculate only sometimes to make it faster
                [fourier_error(iter,g_ind)] = get_fourier_error(modF, aPsi, noise,mask, par.likelihood);
            end
            

        
            if (par.multilayer_object || par.multilayer_probe)
                [Psi(end),R] = modulus_constraint(modF,aPsi,Psi(end), mask, noise, par, 1); %  apply only on the last layer !!! 
            else
                [Psi,R] = modulus_constraint(modF,aPsi,Psi, mask, noise, par,1);
            end
            if iter == 0
                % in the first iteration only find optimal scale for the probe
                probe_amp_corr(1) = probe_amp_corr(1) + Ggather(sum(modF(:).^2));
                probe_amp_corr(2) = probe_amp_corr(2) + Ggather(sum(aPsi(:).^2));
                continue
            end
            if strcmp(par.likelihood, 'poisson')    % calculate only for the first mode
                %% automatically find  optimal step-size, note that for Gauss it is 1 !! 
                cache.beta_xi(g_ind)  =  gradient_descent_chi_solver(self,modF, aPsi2, R,mask, g_ind, mean(cache.beta_xi), cache);
            end
            
            if(par.multilayer_object || par.multilayer_probe)
                ind_modes = par.Nlayers:-1:1;
            else
                ind_modes = 1:max(par.probe_modes, par.object_modes);
            end
            

            for ll = ind_modes
                layer = ll; 
                chi = back_fourier_proj(Psi{min(end,ll)}, self.modes{min(end,ll)})-psi{min(end,ll)};
                               
                 
                %% get optimal gradient lenghts 
                 object_update=0; probe_update=0;m_probe_update= 0;
                
                 if iter >= par.object_change_start && (ll <= max(par.Nlayers, par.object_modes) || par.apply_multimodal_update)
                     object_update = Gfun(grad,chi, probe{min(ll,end)},...
                        probe_max{min(end,ll)}(1,1,min(end,p_ind{ll})),par.delta_p);
                 end
                
                
                if iter >= par.probe_change_start && ll <= max(par.probe_modes, par.Nlayers)
                    %% find optimal probe update !!! 
                        probe_update = Gfun(grad,chi,obj_proj{min(end,ll)},...
                            object_max{min(ll,end)}, par.delta_p);  
                       m_probe_update = mean(probe_update,3);
                end

                if ((ll == 1 && ~(par.multilayer_object || par.multilayer_probe)) || ...
                    (ll == par.Nlayers && (par.multilayer_object || par.multilayer_probe)))  && ...
                     (par.beta_LSQ ||  iter >= par.probe_position_search) 
                    %% variable step extension, apply only the first mode except the 3PIE case 
                    %% it will use the same alpha for the higher modes !!!!    
                                        
                    [cache.beta_probe(g_ind),cache.beta_object(g_ind)] =  gradient_projection_solver(self,chi,obj_proj{min(end,ll)},probe{ll},...
                        object_update, m_probe_update,p_ind{ll}, par, cache);
                     if any(g_ind ==1)
                         verbose(1,'Average alpha p:%3.3g  o:%3.3g ', mean(cache.beta_probe),mean(cache.beta_object));
                     end
                end
                     

                
                
                %%%%%%%%%%%%%%%%%%%%% PROBE UPDATE   %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%
                %% update probe first 
                if ((iter >= par.probe_change_start) && ll <= max(par.probe_modes,par.Nlayers) ...
                        &&  is_method(par, 'PIE')) || ...
                        (par.multilayer_object && ll > 1 && ll <= par.probe_modes)
                        % only in case of first layer probe, otherwise update interprobes
                    
                    beta_probe = get_vals(cache.beta_probe,g_ind) .* get_vals(cache.beta_xi,g_ind);                       
                            
                    if is_method(par, {'ePIE', 'hPIE'})
                        %%%%%%%%%%% update probe %%%%%%%%%%%%%%%%%%%%%%%%%%
                        probe{ll} = Gfun(@upd_probe_Gfun,probe{ll},probe_update, beta_probe); 

                        if (ll == 1 && par.apply_subpix_shift) 
                            probe{ll}  = apply_subpx_shift(probe{ll} , -self.modes{min(end,ll)}.sub_px_shift(g_ind,:));
                        end     
        
                     
                        if iter >= par.probe_change_start  
                            if (par.variable_probe && ll == 1)
                                self.probe{ll}(:,:,p_ind{ll}) = probe{ll};  % slowest line for large datasets !!!!!
                            elseif (~par.multilayer_object || ll == 1) && ll <= par.probe_modes
                                self.probe{ll} = mean(probe{ll},3); % merge information from all the shifted  probes if needed
                            end
                        end
                    end
                end
                             
                if iter > par.probe_fourier_shift_search && ll == 1
                    % search position corrections in the Fourier space, use
                    % only informatiom from the first mode, has to be after
                    % the probes updated 
                    self.modes{1} = gradient_fourier_position_solver(chi, obj_proj{1},probe{1},self.modes{1}, g_ind);
                end
                
                if  iter >= par.probe_position_search
                    % find optimal position shift that minimize chi{1} in current iteration 
                    %[pos_update, cache] = gradient_position_solver(self, chi, obj_proj{1},probe{1,layer}, g_ind, iter, cache);
                    %modified by YJ to prevent bug
                    [pos_update, ~,~,cache] = gradient_position_solver(self, chi, obj_proj{1},probe{1,layer}, g_ind, iter, cache, par);
                    
                    self.modes{1}.probe_positions(g_ind,:)=self.modes{1}.probe_positions(g_ind,:)+pos_update;
                end
                %%%%%%%%%%%%%%%%%%%%% OBJECT UPDATE   %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%
                if iter >= min([par.object_change_start]) && ....
                        ( ll <= max(par.Nlayers, par.object_modes) || par.apply_multimodal_update ) 
                    if ll ~= 1 && ~(par.multilayer_object || par.multilayer_probe) ; continue; end

                    if iter >= par.object_change_start  % && ~(par.multilayer_probe && ll > 1)  % the objects are just empty 
                        beta_object =  get_vals(cache.beta_object,g_ind) .* get_vals(cache.beta_xi,g_ind);
                        if par.share_object
                            obj_ids = 1;  % update only the first object 
                        else
                            obj_ids = unique(scan_ids{ind_ii});  % update only the objects processed in this block 
                        end   
                        
                        if any(beta_object ~= 1)
                            object_update = bsxfun(@times, object_update, beta_object);
                        end
                                     

                        if  is_method(par, 'ePIE')  % use always in nearfield 
                            % classical ePIE, faster constraint application, but it will fail with too high grouping         
                            self.object = set_views(self.object, object_update,layer, obj_ids, g_ind, cache,  scan_ids{ind_ii},[]);
                        elseif is_method(par, 'hPIE')  %% hybrid PIE
                            if par.Nscans == 1 || par.share_object
                                ind_tmp = 1;
                            else
                                ind_tmp = 1+par.object_modes* ((1:par.Nscans)-1);
                            end
                            for kk = ind_tmp
                                obj_illum_sum{kk,layer}(:) = 0;
                                object_upd_sum{kk,layer}(:) = 0;
                            end
                            object_update = bsxfun(@times,  object_update,  aprobe2); % make is more like dumped LSQ solution 
                            [object_upd_sum,obj_illum_sum] = set_views_rc(object_upd_sum,obj_illum_sum, object_update,aprobe2,layer,obj_ids, g_ind, cache, scan_ids{ind_ii},[]);
                            for kk = ind_tmp
                                 self.object{kk,layer} = Gfun(@object_update_Gfun, self.object{kk,layer},object_upd_sum{kk,layer}, obj_illum_sum{kk,layer}, cache.MAX_ILLUM(min(kk,end)));
                            end
                        else
                            error('Unimplemented method %s ', par.method)
                        end
                    end                   
                end
                if ll > 1 && par.multilayer_object
%                     % apply rescaling to make scaling correction 
                    Psi{ll-1} =  probe{ll};
                end
                if ll > 1 && par.multilayer_probe
                    Psi{ll-1} =  obj_proj{ll} +  object_update;
                end  
            end
            
             if check_avail_memory < 0.2 || ~par.keep_on_gpu
                % slow step that is not needed if there is enough memory ,
                % it can slow down almost twice !!!
                 clear Psi R aPsi2 psi probe_update object_update chi 
                 if par.keep_on_gpu
                    warning('Low GPU memory')
                 end
             end
            
        end
       
        

       
       if par.multilayer_object
           self.probe = self.probe(1); 
       end
 
        if par.variable_probe && iter >= par.probe_change_start
            %disp(size(self.probe{1}))
            [self.probe{1}, self.probe_evolution] = apply_SVD_filter(self.probe{1}, par.variable_probe_modes+1,  self.modes{1});
            %disp(size(self.probe{1}))

        elseif par.variable_probe && iter < par.probe_change_start
            self.probe{1} = probe_0; 
        elseif ( ~isempty(self.probe_support)) &&  iter >= par.probe_change_start
            self.probe{1} = apply_probe_contraints(self.probe{1}, self.modes{1});
        end

       if iter == 0
           % apply initial correction for the probe intensity and return
           % it seems to be safer to underestimate the probe amplitude for variable probe method 
           probe_amp_corr = 0.5*sqrt(probe_amp_corr(1) / probe_amp_corr(2)); %% calculate ratio between modF^2 and aPsi^2
           
           %modified by YJ to fix bug with multilayer
           if par.multilayer_object
               self.probe{1} = self.probe{1}*probe_amp_corr;
           else
               
               for ii = 1:par.probe_modes
                   self.probe{ii} = self.probe{ii}*probe_amp_corr;
               end
           end
           verbose(2,'Probe amplitude corrected by %.3g',probe_amp_corr)
           return
       end
       
end

function probe = upd_probe_Gfun(probe,probe_update, alpha_p)
    probe =  probe + alpha_p.*probe_update;
end

function grad = get_grad_flat(chi,proj,max2, delta )
%   ePIE method
    grad =  (1/max2) * chi .* conj(proj) ;
end

function grad = get_grad_lsq(chi,proj,max2, delta )
    % dumped LSQ method  
    aproj = abs(proj) ;
    grad =  chi .*  conj(proj) .* aproj ./ (aproj.^2 + delta*max2)./sqrt(max2);
end

function object = object_update_Gfun(object,object_upd_sum, obj_illum_sum, max)
    object = object +  object_upd_sum ./ (obj_illum_sum+1e-9* max);
end

function array = get_vals(array, ind)
    if isscalar(array)
        return
    else
        array = reshape(array(ind),1,1,[]);
    end
end
