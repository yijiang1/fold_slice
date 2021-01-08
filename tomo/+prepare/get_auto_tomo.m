% [delta_stack_prealign, obj_interf_pos_x, obj_interf_pos_y ] = get_auto_tomo(param_autotomo,surface_calib_file, omnyposfile)
%
% Description:
%
% The function (1) loads the omnyposfile file and determine the scanning
% positions if get_auto_calibration or auto_alignment is 1 and 
% (2) loads the surface_calib_file to give an initial guess for
% the alignemnt array (deltastack) if auto_alignment is 1
%
% Input: 
%
% par. auto_alignment: 0 or 1 (default)
% par. get_auto_calibration: 0 or 1 (default)
% surface_calib_file (mandatory if auto_alignment=1) 
% omnyposfile (mandatory if auto_alignment=1 or get_auto_calibration=1) 
%
% Output:
%
% delta_stack_prealign: used as initial guess for the alignment
% obj_interf_pos_x and obj_interf_pos_y: Object maximum position based on interferometry
%


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


function [delta_stack_prealign, obj_interf_pos_x, obj_interf_pos_y ] = ...
    get_auto_tomo(par,surface_calib_file, omnyposfile, theta, scanstomo)
import beamline.read_omny_pos

import utils.*
import ptycho.*
import beamline.*


obj_interf_pos_x = [];
obj_interf_pos_y = [];
flag_plot = 1;
delta_stack_prealign = [];


%%% To improve: Shifts of the probe are not yet considered here, see
%%% /cSAXS_sxdm_2013_06_omny/matlab/tomo/autotomo_calibration_porous_S00506_S00930.m
if par.auto_alignment ||par.get_auto_calibration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Determine position of first pixel in the reconstructions %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Loading omny_pos for autoalignment')
    for ii = 1:max(size(scanstomo))
        progressbar(ii, max(size(scanstomo)))
        out_orch = read_omny_pos(sprintf(omnyposfile,scanstomo(ii)));
        positions_real = [out_orch.Average_y_st_fzp*1e-6 out_orch.Average_x_st_fzp*1e-6];
%         clear positions
%         positions = positions_real./par.pixel_size;
%         
%         % Change from object to probe positions
%         positions = -positions;
%         
%         positions(:,1) = positions(:,1) - min(positions(:,1));
%         positions(:,2) = positions(:,2) - min(positions(:,2));
%         positions = round(positions);
        
        %%% Object maximum position based on interferometry - sample motion
        %%% Corresponds to pos to coordinates of (1,1) pixel
        %%% increasing number means the sample was higher
        obj_interf_pos_y(ii) = max(positions_real(:,1));
        obj_interf_pos_x(ii) = max(positions_real(:,2));
    end
end

if par.auto_alignment && exist(surface_calib_file, 'file')

    %%% Read calibration file and interpolate correction to these angles
    pos_cal = load(surface_calib_file);
    delta_stack_corr_y_filt = spline(pos_cal.thetasort,pos_cal.delta_stack_corr_y_filt,theta);
    delta_stack_corr_x_filt = spline(pos_cal.thetasort,pos_cal.delta_stack_corr_x_filt,theta);

    %%% Interferometer alignment with mirror surface corrections
    delta_stack_prealign(1,:) = delta_stack_corr_y_filt+obj_interf_pos_y;
    delta_stack_prealign(2,:) = delta_stack_corr_x_filt+obj_interf_pos_x;

    %%% Remove constant term from y alignment
    delta_stack_prealign(1,:) = delta_stack_prealign(1,:)-mean(delta_stack_prealign(1,:));

    %%% Remove sin term from correction in x
    [~,indsort] = sort(theta);
    auxfunc = delta_stack_prealign(2,indsort);
    auxfunc = [auxfunc -auxfunc+auxfunc(end)+auxfunc(1)];
    auxfuncft = fft(auxfunc);
    auxfuncft(3:end-1) = 0;
    auxfunc2 = ifft(auxfuncft);    
    auxfunc3 = auxfunc2(1:end/2);
    delta_stack_prealign(2,indsort) = delta_stack_prealign(2,indsort) - auxfunc3;

    delta_stack_prealign = delta_stack_prealign/par.pixel_size;

    if flag_plot
        figure(1);
        clf;
        subplot(2,1,1)
        plot(theta,obj_interf_pos_y,'.')
        title('Interferometer y position [microns]')
        subplot(2,1,2)
        plot(theta,obj_interf_pos_x,'.')
        title('Interferometer x position [microns]')

        figure(2);
        clf;
        subplot(2,1,1)
        plot(theta,delta_stack_prealign(1,:),'.')
        title('Correction in y [pixels]')
        subplot(2,1,2)
        plot(theta,delta_stack_prealign(2,:),'.')
        title('Correction in x [pixels]')
    end
elseif ~exist(surface_calib_file, 'file')
    warning('Missing surface calibration file %s', surface_calib_file)
end
end
