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
    %{
    %% added by YJ
    %% ADVANCED FLY SCAN - LINE SCAN 
    % interpolate the other modes into new positions 
    pos = self.modes{1}.probe_positions;

    %pos = self.probe_positions_0; %% modified by YJ, use initial position
   
    for ll = 1:par.Nmodes            
        ratio = par.flyscan_dutycycle*(ll-1)/par.Nmodes;
        self.modes{ll}.probe_positions = pos(min((1:self.Npos)+1,self.Npos),:)*ratio + (1-ratio)*pos;
        if ~isempty(jumps)
            % expected step continuation
            self.modes{ll}.probe_positions(jumps,:) = bsxfun(@plus, self.modes{ll}.probe_positions(jumps-1,:),step);
        end
    end
    %}
    
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
        for ll = 1:par.Nmodes
            ratio = par.flyscan_dutycycle*(ll-1)/par.Nmodes;
            self.modes{ll}.probe_positions = pos(min((1:self.Npos)+1,self.Npos),:)*ratio + (1-ratio)*pos;
            if ~isempty(jumps)
                % expected step continuation
                self.modes{ll}.probe_positions(jumps,:) = bsxfun(@plus, self.modes{ll}.probe_positions(jumps-1,:),step);
            end
        end
    end
    
end