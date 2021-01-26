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
    
    switch par.flyscan_trajectory
        case 'line' % line scan w. big jumps
            jumps = diff(self.probe_positions_0);
            %step = median(jumps,1);
            jumps = sum(abs(jumps),2);
            % empirical condition
            jumps = find(jumps > 10*median(jumps));
            %% ADVANCED FLY SCAN - LINE SCAN. Modified by YJ for more general use
            % interpolate the other modes into new positions 
            pos = self.modes{1}.probe_positions;
            %pos(:,1): horizontal positions in pixels
            %pos(:,2): vertical positions in pixels 
            N_lines = length(jumps)+1;
            jumps = [0; jumps; length(pos(:,1))];
            for ll = 1:par.Nmodes
                ratio = par.flyscan_dutycycle*(ll-1)/par.Nmodes;
                pos_temp = pos;
                p_lb = 1;
                for i=1:N_lines
                    x = 1:(jumps(i+1)-jumps(i));
                    x_interp = (x(1)+ratio):1:(x(end)+ratio);
                    p_ub = p_lb + length(x)-1;

                    pos_temp(p_lb:p_ub,1) = interp1(x,pos(p_lb:p_ub,1),x_interp,'pchip');
                    pos_temp(p_lb:p_ub,2) = interp1(x,pos(p_lb:p_ub,2),x_interp,'pchip');  
                    self.modes{ll}.probe_positions = pos_temp;
                    p_lb = p_ub+1;
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
        case 'continuous' % assume that smooth path is used 
            for ii = 1:par.Nscans
                ind = self.reconstruct_ind{ii};
                pos = self.probe_positions_0(ind,:);
                pos_temp = pos;

                for ll = 1:par.Nmodes
                    ratio = par.flyscan_dutycycle*(ll-1)/par.Nmodes;
                    x = 1:size(pos,1);
                    x_interp = (x(1)+ratio):1:(x(end)+ratio);

                    pos_temp(:,1) = interp1(x,pos(:,1),x_interp,'pchip');
                    pos_temp(:,2) = interp1(x,pos(:,2),x_interp,'pchip');  
                    self.modes{ll}.probe_positions(ind,:) = pos_temp;                    
                end
                %{
                % MO's code - seems only work for spiral trajectory
                assert(~any(isfinite(par.probe_position_search)), 'Position refinement and fly scans not suported')

                ind = self.reconstruct_ind{ii};
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
                %}
            end
        case 'external' % use positions from external measurements
    
    
    end
end