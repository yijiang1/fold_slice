% A script to analyze vertica and horizontal through focus scans to 
% determine the size and position of the horizontal and vertical focii

addpath ..
close all
clear
%% vertical beam

scans = [797:806];

p.motor_name = 'py';
p.title_str = 'Vertical beam';
p.motor_units = '(microns)';
p.counter = 'diode';
p.pausetime = 0.5;
p.coarse_motor = 'hz';

out_ver = utils.focus_series_fit(scans,p);

%% horizontal beam

scans = [650:660];%[183:193];%[505:514];% scans 108 to, 89 to 

p.motor_name = 'px';
p.title_str = 'Horizontal beam';
p.motor_units = '(microns)';
p.counter = 'diode';
p.pausetime = 0.5;
p.coarse_motor = 'hz';

out_hor = utils.focus_series_fit(scans,p);

%% plot both
figure(7)
plot(out_ver.fitout,'b',out_ver.coarse_motor,out_ver.fwhm,'bo');
hold on
plot(out_hor.fitout,'r',out_hor.coarse_motor,out_hor.fwhm,'ro');
title(sprintf('beam focus, vertex (hor,ver) (%.1f, %.1f)',out_hor.vertex,out_ver.vertex))
xlabel(p.coarse_motor)
ylabel(sprintf('FWHM %s',p.motor_units))
legend('vertical','fit','horizontal','fit')
hold off
fprintf('The vertical vertex of the parabola is at %s = %f\n',p.coarse_motor,  out_ver.vertex);
fprintf('The horizontal vertex of the parabola is at %s = %f\n',p.coarse_motor,out_hor.vertex);

%% Example for sgalil continuous scan

scans = [797:806];

p.position_file = '~/Data10/sgalil/S%05d.dat';
p.fast_axis_index = 1;  % = 1 or 2 for x and y scan respectively
p.mcs_file = sprintf('~/Data10/mcs/S%02d000-%02d999/S%%05d/%s_%%05d.dat',floor(scans(1)/1000),floor(scans(1)/1000),beamline.identify_eaccount);
p.mcs_channel = 3;
p.title_str = 'Horizontal beam';
p.motor_units = '(mm)';
p.counter = 'diode';
p.pausetime = 0.5;
p.coarse_motor = 'samy';

out_ver = utils.focus_series_fit(scans,p);

%%

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2018 by Paul Scherrer Institute (http://www.psi.ch)    |
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