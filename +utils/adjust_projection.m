% function used to correct the image orientation of mcs_mesh data.
% [output, output_pos] = adjust_projection(input, snake_scan, fast_axis_x, positions)
%       input = data to be corrected. For mcs the data should be a 2D matrix.
%       snake_scan = 0 for off and 1 for on
%       fast_axis_x = 1 for fast axis along x, 0 for fast axis along y
%
%       output = corrected data
%       output_pos = corrected output positions, could be used to see if
%           there was a problem with the correction
% For snake scans the routine decides the flipping based on the positions
 
%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
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

function [output, output_pos] = adjust_projection(input, snake_scan, fast_axis_x, positions)

if nargin < 4
    positions = [];
end

output_temp = input;
output_pos  = positions;

%%% Sanity checks %%%
if ~isempty(output_pos)
    %%% fast axis direction %%%
    average_x_step_fast_axis = mean(mean(abs(diff(output_pos(:,:,1),1,1))));
    average_y_step_fast_axis = mean(mean(abs(diff(output_pos(:,:,2),1,1))));
    fast_axis_x_from_pos = fast_axis_x;
    if (fast_axis_x)&&(average_x_step_fast_axis < average_y_step_fast_axis)
        warning('You specified fast_axis_x true, but the positions seem to be for fast axis along y')
        fast_axis_x_from_pos = false;
        fast_axis_ind = 2;
        slow_axis_ind = 1;
    elseif (~fast_axis_x)&&(average_x_step_fast_axis > average_y_step_fast_axis)
        warning('You specified fast_axis_x false, but the positions seem to be for fast axis along x')
        fast_axis_x_from_pos = true;
        fast_axis_ind = 1;
        slow_axis_ind = 2;
    elseif fast_axis_x
        fast_axis_ind = 1;
        slow_axis_ind = 2;
    elseif ~fast_axis_x
        fast_axis_ind =2;
        slow_axis_ind = 1;
    end
    %%% Scan quality check %%%
    if ~any(output_pos(:)==0)
        aux_fast = abs(diff(output_pos(:,:,fast_axis_ind),1,1));
        aux_slow = abs(diff(output_pos(:,:,slow_axis_ind),1,2));
        average_fastaxis_absstep = mean(aux_fast(:));
        average_slowaxis_absstep = mean(aux_slow(:));
        std_fastaxis_absstep = std(aux_fast(:));
        std_slowaxis_absstep = std(aux_slow(:));
        step_text = sprintf('\n Step, (fast axis,slow axis) +/- (std,std) = (%.2f,%.2f) +/- (%.2f,%.2f) microns.',...
            average_fastaxis_absstep*1e3,average_slowaxis_absstep*1e3,std_fastaxis_absstep*1e3,std_slowaxis_absstep*1e3);
        if (std_fastaxis_absstep>average_fastaxis_absstep*0.05)||(std_slowaxis_absstep>average_slowaxis_absstep*0.05)
            warning(step_text)
            pause(2)
        else
            if nargout>1
                disp(step_text);
            end
        end
    end
    
    

    
    %%% snake scans %%%
    average_fastaxis_step = mean(diff(output_pos(:,:,fast_axis_ind),1,1));
    % is this a snake scan?
    if abs(average_fastaxis_step(2)-average_fastaxis_step(3))==0
        % Do nothing, data has not been loaded
        snake_scan_from_pos = snake_scan;
    elseif abs(average_fastaxis_step(2)-average_fastaxis_step(3))>abs(average_fastaxis_step(1))
        snake_scan_from_pos = true;
    else
        snake_scan_from_pos = false;
    end
    
    if snake_scan ~= snake_scan_from_pos
        warning(['You specified snake_scan = ' num2str(snake_scan) ' but from the positions it seems that snake_scan = ' num2str(snake_scan_from_pos)])
    end
end



%%% Handling the flipping of the data %%%
if snake_scan %
    startflipind = 1;
    if ~isempty(output_pos)
        if mean(diff(output_pos(:,1,fast_axis_ind),1,1))>0
            startflipind = 2;
        else
            startflipind = 1;
        end
        output_pos(:,startflipind:2:end,:)  = flipud(output_pos(:,startflipind:2:end,:));
    end
    output_temp(:,startflipind:2:end,:,:) = flipud(output_temp(:,startflipind:2:end,:,:));
end

if fast_axis_x
    output_temp = permute(output_temp,[2 1 3]);
    output_pos  = permute(output_pos, [2 1 3]);
end

output     = rot90(output_temp,2);
output_pos = rot90(output_pos,2);


% if ~isempty(output_pos)
%     if (mean(mean(diff(output_pos(:,:,1),1,2)))>0)||(mean(mean(diff(output_pos(:,:,2),1,2)))>0)
%         warning('Something is wrong with the positions, I dont know what so Ill show you in a figure of the positions after adjusting them. Positions should monotonically decrease with increase x or y coordinate')
%         figure(123)
%         subplot(1,2,1)
%         imagesc(output_pos(:,:,1))
%         title('X position')
%         subplot(1,2,2)
%         imagesc(output_pos(:,:,2))
%         title('Y position')
%     end
% end
    
    
return    
    

