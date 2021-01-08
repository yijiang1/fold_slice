% GRADIENT_POSITION_SOLVER solve position errors in the real space 
%
% [pos_update, cache] = gradient_position_solver(self,xi,O,P,ind, iter, cache)
%
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** xi        exit wave update vector 
% ** O         object views 
% ** P         probe or probes 
% ** ind       indices of of processed position 
% ++ iter      current iteration 
% ** cache     structure with precalculated values to avoid unnecessary overhead
%
% returns:
% ++ pos_update     position updates for each of the indices   
% ++ cache          updated  structure with precalculated values 
%
% see also: engines.GPU_MS.LSQML, engines.GPU_MS.PIE 


function [pos_update, probe_rotation,probe_scale,cache] = gradient_position_solver(self,xi,O,P,ind, iter, cache, par)
    import engines.GPU_MS.GPU_wrapper.*
    import math.*
    import utils.*
    % use gradinent solver for position correction  

    low_mem_errs = {'id:parallel:gpu:array:OOMForOperation',...
     'id:MATLAB:LowGPUMem','MATLAB:LowGPUMem',...
     'parallel:gpu:array:OOM',...
     'parallel:gpu:device:UnknownCUDAError', ...
     'parallel:gpu:array:OOMForOperation', ...
     'parallel:gpu:array:FFTInternalError'};
    
    % wrapper around get_img_grad, in case of low memory it will try to repeat
    % Ntimes before giving up 
    pos_update = 0; probe_rotation = 0; probe_scale = 0; 
    N = 5; 
    for ii = 1:N
        try
            % reuse dx_O, dy_O to save memory !! 
            [dx_O,dy_O]=get_img_grad(O);
            
            if iter >= par.detector_rotation_search
                %% estimate detector rotation 
                xgrid  = Garray(linspace(-1,1,self.Np_p(1))'); 
                ygrid  = Garray(-linspace(-1,1,self.Np_p(2))); 
                [nom, denom] = Gfun(@get_coefs_mixed,xi, P, dx_O, dy_O, xgrid, ygrid);
                probe_rotation =  gather(sum2(nom)./ sum2(denom));
            end
            
            if iter >= par.detector_scale_search
                %% estimate detector scale (ie pixel scale error in farfield mode) 
                xgrid  = Garray(-linspace(-1,1,self.Np_p(2)) .* tukeywin(self.Np_p(2), 0.1)'); 
                ygrid  = Garray(-linspace(-1,1,self.Np_p(1))'.* tukeywin(self.Np_p(1), 0.1)); 
                [nom, denom] = Gfun(@get_coefs_mixed,xi, P, dx_O, dy_O, xgrid, ygrid);
                probe_scale =  gather(sum2(nom)./ sum2(denom));
                probe_scale = 0.5*mean(probe_scale) / mean(self.Np_p);
            end

            if iter >= par.probe_position_search
                %% estimate sample shift 
                [dx_O, denom_dx, dy_O, denom_dy] = Gfun(@get_coefs_shift,xi,P,dx_O, dy_O);
                dx =  sum2(dx_O)./ sum2(denom_dx);
                dy =  sum2(dy_O)./ sum2(denom_dy);
            end
            break
        catch ME
            warning('Low memory')
            if ~any(strcmpi(ME.identifier, low_mem_errs))
                rethrow(ME)
            end
            pause(1)
        end
    end
    if ii == 5
        rethrow(ME) 
    end

    if iter < par.probe_position_search
        return
    end
    
    shift = squeeze(Ggather(cat(4,dx, dy)));
    
    % modified by YJ. allow user to spcify maximum position update
    max_shift = min(par.max_pos_update_shift, 10*mad(shift)); %why a factor of 10??
    
    % prevent outliers and too rapid shifts 
    %max_shift = min(0.1, 10*mad(shift)); 
    shift = min(abs(shift), max_shift) .* sign(shift); % avoid too fast jumps, <0.5px/iter is enough 
    %old code
    %shift = min(abs(shift), 0.2) .* sign(shift); % avoid too fast jumps, <0.5px/iter is enough 

    pos_update = reshape(shift,[],2); 
    
    if ~isfield(cache, 'velocity_map_positions')
        cache.velocity_map_positions = zeros(self.Npos,2,'single');
    end
    if ~isfield(cache, 'position_update_memory')
        cache.position_update_memory = {};
    end
    
    cache.position_update_memory{iter}(ind,:) = pos_update;
   
    %% USE MOMENTUM ACCELERATION TO MAKE THE CONVERGENCE FASTER
    ACC = 0; 
    if par.probe_position_search_momentum >0
    	momentum_memory = par.probe_position_search_momentum; % 
     
         % only in case far field ptychography 
        if isinf(self.z_distance) && sum(cellfun(@length, cache.position_update_memory) > 0) > momentum_memory 
            %corr_level = zeros(momentum_memory,1);
            for ii = 1:momentum_memory
            	corr_level(ii) = mean(diag(corr(cache.position_update_memory{end}(ind,:), cache.position_update_memory{end-ii}(ind,:))));
            end
            if all(corr_level > 0 )
            	%estimate optimal friction from previous steps 
	        	poly_fit = polyfit(0:momentum_memory,log([1,corr_level]),1); 
 
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 gain = 0.5;                           % smaller -> lower relative speed (less momentum)
                 friction =  0.1*max(-poly_fit(1),0);   % smaller -> longer memory, more momentum 
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             else
                gain = 0; friction = 0.5; 
            end
         
            cache.velocity_map_positions(ind,:) = cache.velocity_map_positions(ind,:)*(1-friction) + pos_update;
	         % apply the velocity to the refined positions , if the postition updated are sufficiently small 

             if max(abs(pos_update)) < 0.1
                 ACC = norm2(pos_update + gain*cache.velocity_map_positions(ind,:)) / norm2(pos_update); 
                 pos_update = pos_update + gain*cache.velocity_map_positions(ind,:); 
             end 
        end
    end
%     try
    %ACC = 0; 
%     momentum_memory = 5; % remember 5 iterations 
%     
%     % only in case far field ptychography 
%     if isinf(self.z_distance) && sum(cellfun(@length, cache.position_update_memory) > 0) > momentum_memory 
%         for ii = 1:momentum_memory
%             corr_level(ii) = mean(diag(corr(cache.position_update_memory{end}(ind,:), cache.position_update_memory{end-ii}(ind,:))));
%         end
%         if all(corr_level > 0 )
%             %estimate optimal friction from previous steps 
%             poly_fit = polyfit(0:momentum_memory,log([1,corr_level]),1); 
% 
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             gain = 0.5;                           % smaller -> lower relative speed (less momentum)
%             friction =  0.1*max(-poly_fit(1),0);   % smaller -> longer memory, more momentum 
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         else
%            gain = 0; friction = 0.5; 
%         end
%         
%         cache.velocity_map_positions(ind,:) = cache.velocity_map_positions(ind,:)*(1-friction) + pos_update;
%         % apply the velocity to the refined positions , if the postition updated are sufficiently small 
% 
%         if max(abs(pos_update)) < 0.1
%             ACC = norm2(pos_update + gain*cache.velocity_map_positions(ind,:)) / norm2(pos_update); 
%             pos_update = pos_update + gain*cache.velocity_map_positions(ind,:); 
%         end
% 
%     end
%     catch
%         keyboard
%     end
    if any(ind==1)
        verbose(1,'Grad pos corr -- AVG step  %3.3g px , acceleration = %4.1f', max(abs(pos_update(:))), ACC)
    end
    
end
    
function [nom1, denom1, nom2, denom2] = get_coefs_shift(xi, P, dx_O, dy_O)

    dx_OP = dx_O.*P;
    nom1 = real(conj(dx_OP) .* xi);
    denom1 = abs(dx_OP).^2; 

    dy_OP = dy_O.*P;
    nom2 = real(conj(dy_OP) .* xi);
    denom2 = abs(dy_OP).^2; 

end

    
function [nom, denom] = get_coefs_mixed(xi, P, dx_O, dy_O, xgrid, ygrid)

    dm_O = dx_O .* xgrid + dy_O .* ygrid; 

    dm_OP = dm_O.*P;
    nom = real(conj(dm_OP) .* xi);
    denom = abs(dm_OP).^2; 

end
