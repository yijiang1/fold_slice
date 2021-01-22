% PREPARE_FLYSCAN_POSITIONS from finit number of measured position interpolate possitions for each
% measured frame when fly scan is used
% 
% self = prepare_flyscan_positions(self, par)
% 
% 
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** par       structure containing parameters for the engines 
%
% returns: 
% ++ self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..


function self = prepare_flyscan_positions(self, par)
import engines.GPU.GPU_wrapper.*
import math.*
import utils.*
import plotting.*

    jumps = diff(self.probe_positions_0);
    step = median(jumps,1);
    jumps = sum(abs(jumps),2);
    % empirical condition 
    jumps = find(jumps > 10*median(jumps));
    
    if length(jumps) < par.Nscans
        % assume that smooth path is used 
       %% ADVANCED FLY SCAN - SPIRAL
        for ii = 1:par.Nscans
            assert(~any(isfinite(par.probe_position_search)), 'Position refinement and fly scans not suported')

            ind =  self.reconstruct_ind{ii};
            pos = self.probe_positions_0(ind,:);
            [ang, rad] = cart2pol(pos(:,1)-pos(1,1), pos(:,2)-pos(1,2)); 
            ang = unwrap(ang); 
            % get interpolate d positions of the sub probes 
            ang_all = ang + (par.flyscan_offset -0.5+linspace(0,par.flyscan_dutycycle*(par.Nmodes-1)/par.Nmodes, par.Nmodes)   ).*[diff(ang);0]; 
            rad_all = interp1(ang, rad, ang_all, 'pchip'); 
            [X,Y] = pol2cart(ang_all, rad_all);
            for ll = 1:par.Nmodes
                    self.modes{ll}.probe_positions(ind,:) = pos(1,1:2) + [X(:,ll), Y(:,ll)];
                %if iter == 1;  self.probe{ll} = self.probe{1}; end
            end
        end
    else
        %% ADVANCED FLY SCAN - LINE SCAN 
        %disp('FLY SCAN - LINE SCAN')
        % interpolate the other modes into new positions 
        pos = self.modes{1}.probe_positions;
        %modified by YJ
        %pos(:,1): horizontal positions in pixels
        %pos(:,2): vertical positions in pixels
        %TODO: suport non-rectangular grid
        for ll = 1:par.Nmodes
            ratio = par.flyscan_dutycycle*(ll-1)/par.Nmodes;
            pos_temp = pos;
            Ny = length(jumps)+1;
            Nx = length(pos)/Ny;
            x = 1:Nx;
            x_interp = (1+ratio):1:(Nx+ratio);
            
            for i=1:Ny
                p_lb = (i-1)*Nx+1;
                p_ub = i*Nx;
                pos_temp(p_lb:p_ub,1) = interp1(x,pos(p_lb:p_ub,1),x_interp,'pchip');
                pos_temp(p_lb:p_ub,2) = interp1(x,pos(p_lb:p_ub,2),x_interp,'pchip');  
                self.modes{ll}.probe_positions = pos_temp;
            end
        end
        
        %{ 
        %MO's code - may have some bugs
        for ll = 1:par.Nmodes
            ratio = par.flyscan_dutycycle*(ll-1)/par.Nmodes;
            self.modes{ll}.probe_positions = pos(min((1:self.Npos)+1,self.Npos),:)*ratio + (1-ratio)*pos;
            if ~isempty(jumps)
                % expected step continuation
                self.modes{ll}.probe_positions(jumps,:) = bsxfun(@plus, self.modes{ll}.probe_positions(jumps-1,:),step);
            end
        end
        %}
    end
    
end