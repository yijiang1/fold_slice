function [ppX_rot,ppY_rot]=scan_position_rot(N_scan_x,N_scan_y,scanStepSize_x,scanStepSize_y,rot_ang)
% scan positions after rotation
    ppx = linspace(-floor(N_scan_x/2),ceil(N_scan_x/2)-1,N_scan_x)*scanStepSize_x;
    ppy = linspace(-floor(N_scan_y/2),ceil(N_scan_y/2)-1,N_scan_y)*scanStepSize_y;
    [ppX,ppY] = meshgrid(ppx,ppy);

    ppY_rot = ppX*(-sind(rot_ang)) + ppY*cosd(rot_ang);
    ppX_rot = ppX*cosd(rot_ang) + ppY*sind(rot_ang);

end
