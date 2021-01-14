% function [ out ] = focus_series_fit( scans, p )
% Receives scan numbers and parameters as a structure p
% Input:
%   scans
%       For SPEC variables
%   p.motor_name    From SPEC
%   p.counter       From SPEC
%       For sgalil position file
%   p.position_file     Example  '~/Data10/sgalil/S%05d.dat'
%   p.fast_axis_index   (= 1 or 2) for x or y scan respectively
%       For mcs counter
%   p.mcs_file          Example sprintf('~/Data10/mcs/S%02d000-%02d999/S%%05d/%s_%%05d.dat',floor(scans(ii)/1000),floor(scans(ii)/1000),beamline.identify_eaccount);
%   p.mcs_channel       Channel number, e.g. = 3
%
% Optional
%   p.motor_units
%   p.plot
%   p.title_str 
%   p.coarse_motor
%
% Output
%   out.fitout        Parameters of quadratic fit
%   out.coarse_motor  Coarse motor name is passed back
%   out.fwhm          A vector with the fwhm for each scan
%   out.vertex        The position of  coarse motor with minimum fwhm from the quadratic fit


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

function [ out ] = focus_series_fit( scans, p )


out = struct;

if isempty(scans)
    error('Scans input seems to be empty')
end
if ~isfield(p,'plot')
    p.plot = true;
end
if ~isfield(p,'title_str')
    p.title_str = '';
end
if ~isfield(p,'motor_units')
    p.motor_units = '';
end
if ~isfield(p,'coarse_motor')
    p.motor_units = '';
end
if ~isfield(p,'pausetime')
    p.pausetime = 0;
end
% mcs
if ~isfield(p,'mcs_file')
    p.mcs_file = [];
end
if ~isfield(p,'mcs_channel')
    p.mcs_channel = [];
end

% sgalil
if ~isfield(p,'position_file')
    p.position_file = [];
end
if ~isfield(p,'fast_axis_index')
    p.fast_axis_index = 1;
end

S_all=io.spec_read('~/Data10/','ScanNr',scans);
width = scans*0;
coarse_motor = scans*0;

for ii=1:length(scans)
    
    if numel(S_all) == 1
        S{1} = S_all;
    else
        S = S_all;
    end
    
    if isempty(p.mcs_file)
        y = getfield(S{ii},p.counter); %#ok<GFLD>
        y(1:end-1)=diff(y);
        y(end) = 0;
        y(end)=y(end-1);
    else
        data = io.image_read(sprintf(p.mcs_file,scans(ii),scans(ii)));
        y = squeeze(data.data(p.mcs_channel,1,:));
        y(1:end-1)=diff(y);
        y([end end+1]) = 0;
    end

    if isempty(p.position_file)
        x = getfield(S{ii},p.motor_name); %#ok<GFLD>
    else
        data = io.image_read(sprintf(p.position_file,scans(ii)));
        x = data.data(p.fast_axis_index,:).';
    end
    
%        General model Gauss1:
%      f(x) =  a1*exp(-((x-b1)/c1)^2)
%      Coefficients (with 95% confidence bounds):
%        a1 =       -2754  (-2839, -2669)
%        b1 =      -84.29  (-84.29, -84.28)
%        c1 =    0.002197  (0.002118, 0.002276)
    
    [yabsmax, ind_absmax] = max(abs(y));
%     p0.a1 = y(ind_absmax);
%     p0.b1 = x(ind_absmax);
%     p0.c1 = 1e-9;
    p0 = [y(ind_absmax) x(ind_absmax) 1e-3];

%     f = fit(x,y,'gauss1');
    f = fit(x,y,'gauss1', 'StartPoint', p0 );

    if p.plot
        figure(4)
        plot(f,x,y,'.-');
        title(p.title_str)
        xlabel(sprintf('%s %s',p.motor_name,p.motor_units))
        ylabel(p.counter)
        drawnow
        pause(p.pausetime)
    end

    width(ii)=f.c1*2*sqrt(2*log(2))/sqrt(2);
    fprintf('S%05d, FWHM = %.2e %s\n',scans(ii),width(ii),p.motor_units)
    coarse_motor(ii)=getfield(S{ii},p.coarse_motor); %#ok<GFLD>
end

figure(5)
plot(coarse_motor,width,'-bo')
title(p.title_str)
xlabel(p.coarse_motor)
ylabel(sprintf('FWHM %s',p.motor_units))

if numel(scans)>2
    
    h = fit(coarse_motor.',width.','poly2');
    
    figure(6)
    plot(h,coarse_motor,width);
    title(p.title_str)
    xlabel(p.coarse_motor)
    ylabel(sprintf('FWHM %s',p.motor_units))
    
    vertex = -h.p2/(2*h.p1);
    fprintf('\n\nThe vertex of the parabola is at %s = %f\n\n',p.coarse_motor,vertex)
    
    fprintf('Average FWHM = %.2e %s\n',mean(width),p.motor_units)
    fprintf('Minimum FWHM = %.2e %s\n',min(width),p.motor_units)
    fprintf('Maximum FWHM = %.2e %s\n',max(width),p.motor_units)
    
    out.fitout = h;
    out.coarse_motor = coarse_motor;
    out.fwhm = width;
    out.vertex = vertex;
end

end

