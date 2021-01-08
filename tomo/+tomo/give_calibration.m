%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: 
%
% [theta_corr_sorted,delta_stack_corr_x_filt,delta_stack_corr_y_filt] = give_calibration(obj_interf_pos_x, obj_interf_pos_y, deltastack, deltaslice, param)
%
% Description:
%
% The function (1) saves the alignment arrays for later use and 
% (2) saves the vertical correction into a .txt (and .mat) file for
% correcting the vertical fluctuarion throguh SPEC
%
% Input: 
%
% obj_interf_pos_x and obj_interf_pos_y: obatined from running get_auto_tomo.m
% deltastack: alignment for x and y
% deltaslice: additional alignment for x
% param. savedata: 0 (default) or 1
% param. surface_calib_file
% param. get_auto_calibration: 0 or 1 (default) 
% param. pixsize (mandatory) 
% param. theta (mandatory) 
% param. scans (mandatory) 
% param. output_folder (default: './') 
%
% Output:
%
% theta_corr_sorted
% delta_stack_corr_x_filt 
% delta_stack_corr_y_filt
%
% 2017-03-30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


function [thetasort,delta_stack_corr_x_filt,delta_stack_corr_y_filt] = ...
    give_calibration(obj_interf_pos_x, obj_interf_pos_y, deltastack, deltaslice, theta, scans, param)



if isfield(param,'get_auto_calibration')
    get_auto_calibration = param.get_auto_calibration;
else
    fprintf('Using get_auto_calibration = 1\n');
    get_auto_calibration = 1;
end 

if ~isfield(param,'surface_calib_file') 
    error('Set param.surface_calib_file');
end 

if isfield(param,'output_folder')  
    output_folder = param.output_folder;
else
    fprintf('Set output_folder to the pwd \n');
    output_folder = './';
end 

if isfield(param,'pixel_size')
    pixsize = param.pixel_size;
else
    error('Need to specify param.pixel_size\n');
end 


[thetasort, indsort] = sort(theta);
if get_auto_calibration
    Nangles = length(theta);
    
    %%% Get correction from interferometer and alignment values in meters
    delta_stack_corr_y = deltastack(1,:)*pixsize;
    delta_stack_corr_x = (deltastack(2,:)+deltaslice)*pixsize;



    
    %% %%% remove outliers by median filter %%%%%%%%%%%%%%%%%%%%
    
    %%% Subtract the a*sin(x+b)+c term from x correction 
    [rigid_shift] = fit_sinus(theta, delta_stack_corr_x');

    delta_stack_corr_x = delta_stack_corr_x - rigid_shift'; 
    %%% Remove constant term from y correction
    delta_stack_corr_y = delta_stack_corr_y - mean(delta_stack_corr_y);
   
    %%% filtering to avoid outliers in the alignment, seems to work better when it is done
    %%% after sinus removal 
    delta_stack_corr_x_filt = medfilt1(delta_stack_corr_x, 5);
    delta_stack_corr_y_filt = medfilt1(delta_stack_corr_y, 5);

    delta_stack_corr_x_filt = delta_stack_corr_x_filt + rigid_shift'; 
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
 
    %%% apply sorting
    delta_stack_corr_x_filt = delta_stack_corr_x_filt(indsort)';
    delta_stack_corr_y_filt = delta_stack_corr_y_filt(indsort)';
    
    figure(1);
    clf;
    subplot(2,1,1)
    plot(theta,obj_interf_pos_y*1e6,'.')
    title('Interferometer y position')
    axis tight ; grid on 
    xlabel('Angles [deg]')
    ylabel('Shift [\mum]')
    subplot(2,1,2)
    plot(theta,obj_interf_pos_x*1e6,'.')
    title('Interferometer x position (plus "clicking correction" term) [microns]')
    axis tight ; grid on 
    xlabel('Angles [deg]')
    ylabel('Shift [\mum]')
        
    figure(4);
    clf;
    subplot(2,1,1)
    plot(theta,deltastack(1,:)*pixsize*1e6,'.')
    title('Alignment in y (without removal of constant term)')
    ylabel('Shift [\mum]')
    xlabel('Angles [deg]')
    axis tight ; grid on 

    subplot(2,1,2)
    plot(theta,( deltastack(2,:)+deltaslice)*pixsize*1e6,'.')
    title('Alignment in x (without removal of sinus term)')
    ylabel('Shift [\mum]')
    xlabel('Angles [deg]')
    axis tight ; grid on 

    figure(3);
    clf;
    subplot(2,1,1)
    plot(thetasort,delta_stack_corr_y_filt*1e6,'.')
    title('Correction in y (after removal of constant term)')
    ylabel('Shift [\mum]')
    axis tight ; grid on 
    xlabel('Angles [deg]')
    subplot(2,1,2)
    plot(thetasort,delta_stack_corr_x_filt*1e6,'.')
    title('Correction in x (after removal of constant term)')
    ylabel('Shift [\mum]')
    axis tight ; grid on 
    xlabel('Angles [deg]')
    output_png1 = [ output_folder 'delta_stack_corr.png'];
    fprintf('Writting image files \n %s\n',output_png1);
    print('-f3','-dpng','-r300',output_png1);

    

   
    %% Calculate the correction for interferometers 
    corr(:,1) =  delta_stack_corr_x_filt*1e6;
    corr(:,2) =  delta_stack_corr_y_filt*1e6;
    
    %%% Filter the correction to remove effects of long term drifts !!
    corr = medfilt1(corr, 11);
    filter = ceil(Nangles / 10); 
    corr(:,1) = smooth(thetasort, corr(:,1), filter, 'sgolay');
    corr(:,2) = smooth(thetasort, corr(:,2), filter, 'sgolay');
    
    calibration_file = sprintf('correction_interferometers_um_S%05d.txt',scans(1)); 


    % -----
    figure(50);
    clf()
    plot(thetasort,corr, '.-'); grid on; ylabel('Correction \mum')
    legend({'Horizontal', 'Vertical'})
    title(sprintf('%s',output_folder),'interpreter','none');
    xlabel('Angles [deg]')
    axis tight; grid on 
    plotting.suptitle(sprintf('Final interferometer correction in file \n%s', calibration_file),'interpreter', 'none')

    output_png1 = fullfile( output_folder, 'y_alignment.png');
    fprintf('Writting image files \n %s\n',output_png1);
    print('-f50','-dpng','-r300',output_png1);


    figure(51);
    plot(1:length(theta),theta); grid on; 
    title(sprintf('%d projections',length(theta)));

    output_png1 = fullfile(  output_folder , 'theta.png');
    fprintf('Writting image files \n %s\n',output_png1);
    print('-f51','-dpng','-r300',output_png1);


    %%% Save file
    utils.verbose(0, 'Saving interferometer correction to %s', param.surface_calib_file)
    obj_interf_pos_x_sort = obj_interf_pos_x(indsort); 
    obj_interf_pos_y_sort = obj_interf_pos_y(indsort); 
    utils.savefast_safe(param.surface_calib_file, 'thetasort','delta_stack_corr_x_filt','delta_stack_corr_y_filt', 'obj_interf_pos_x_sort', 'obj_interf_pos_y_sort')
    
    
    %%% save correction for SPEC 
    if ~exist(calibration_file, 'file') || strcmpi(input(sprintf('Do you want to overwrite %s (y/N)? ', calibration_file),'s'),'y')
        utils.verbose(0, 'Saving interferometer correction to %s', calibration_file)

        h = fopen(calibration_file,'w');
        fprintf(h,'corr_elements = %d \n', Nangles);
        fprintf(h,'corr_elements_x = %d \n', Nangles);
        for jj = 1:Nangles       
            fprintf(h,'%s[%d] = %.6f \n',...
                'corr_angle',jj-1,thetasort(jj));
            fprintf(h,'%s[%d] = %.6f \n',...
                'corr_angle_x',jj-1,thetasort(jj));
            fprintf(h,'%s[%d] = %.6f \n',...
                'corr_pos',jj-1, corr(jj,2));
            fprintf(h,'%s[%d] = %.6f \n',...
                'corr_pos_x',jj-1, corr(jj,1));
        end
        fclose(h);
    
    end
    
end
end


function [rigid_shift] = fit_sinus(theta, signal)
    % subtract the a*sin(x+b) from the data 
     orthbase = [sind(theta(:)), cosd(theta(:)),ones(length(theta),1)]; % 
     coefs = (orthbase'*orthbase) \ (orthbase'*signal);
     % avoid object drifts within the reconstructed FOV
     coefs(3) = 0; % preserve the offset 
     rigid_shift = orthbase*coefs; 
end