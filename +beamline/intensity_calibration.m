% Description:
% Return beam intensity in photons / sec based on calibration with
% a glassy carbon sample and and air
%
% Dependencies:
% spec_read

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
%   using the “cSA% Description:
% Return beam intensity in photons / sec based on calibration with
% a glassy carbon sample and and air
%
% Dependencies:
% spec_readXS matlab package” developed by the CXS group,
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

function [ diodescale ] = intensity_calibration(specfile, air_scanno, gc_scanno, det_dist_mm, varargin)
import io.spec_read
import utils.find_files

pixel_size_mm = 0.172;
gc_file = 'glassycarbon_L14_xsection.dat';
binned_path = '~/Data10/analysis/radial_integration/';

if (nargin < 4)
    fprintf('\nUsage:\n');
    fprintf('%s(specfile, air_scanno, gc_scanno, det_dist_mm, [[,<name>,<value>] ...]);\n\n',mfilename);
    fprintf('specfile is the full path to the SPEC dat-file.\n');
    fprintf('air_scanno and gc_scanno are SPEC scan numbers for empty and Glassy Carbon L14 measurements.\n');
    fprintf('det_dist_mm is the sample to detector distance in mm.\n');
    fprintf('\nThe optional <name>,<value> pairs are:\n');
    fprintf('''PixelSize_mm'', <value in mm>    Size of detector pixel in mm, default is %s\n', pixel_size_mm);
    fprintf('''CrossSectionFile'', <filename>   Full path to file containing the cross section of the standard, default is ''%s''\n', gc_file);
    fprintf('''BinnedPath'', <filepath>         Directory containing the radially binned detector frames, default is ''%s''\n', binned_path);
    fprintf('\n');
    error('Not enough input arguments.');
end

% accept cell array with name/value pairs as well
no_of_in_arg = nargin;
if (nargin == 5)
    if (isempty(varargin))
        % ignore empty cell array
        no_of_in_arg = no_of_in_arg -1;
    else
        if (iscell(varargin{1}))
            % use a filled one given as first and only variable parameter
            varargin = varargin{1};
            no_of_in_arg = 1 + length(varargin);
        end
    end
end

% check number of input arguments
if (rem(no_of_in_arg,2) ~= 0)
    error('The optional parameters have to be specified as ''name'',''value'' pairs');
end

% parse the variable input arguments
vararg = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'PixelSize_mm'
            pixel_size_mm = value;
        case 'CrossSectionFile'
            gc_file = value;
        case 'BinnedPath'
            binned_path= value;
        otherwise
            vararg{end+1} = name; %#ok<AGROW>
            vararg{end+1} = value; %#ok<AGROW>
    end
end


s_air = spec_read(specfile, 'ScanNr', air_scanno);
s_gc = spec_read(specfile, 'ScanNr', gc_scanno);

[dd, air_intfile] = find_files(sprintf(strcat(binned_path, 'e*_1_%05d_00000_00000_integ.mat'), air_scanno));
airint = load(strcat(dd, air_intfile.name));
[dd, gc_intfile]  = find_files(sprintf(strcat(binned_path, 'e*_1_%05d_00000_00000_integ.mat'), gc_scanno));
gcint = load(strcat(dd, gc_intfile.name));

lambda = 12.39852 / s_gc.mokev;
q_gc = 4*pi * sin(0.5*atan(gcint.radius*pixel_size_mm/det_dist_mm)) / lambda;
gc_transmission = (sum(s_gc.diode)/sum(s_gc.sec)) / (sum(s_air.diode)/sum(s_air.sec));
gc_time = sum(s_gc.sec);
I_gc = sum(sum(gcint.I_all, 3), 2) - (sum(s_gc.diode)/sum(s_air.diode)) * sum(sum(airint.I_all, 3), 2);
Ierr_gc = sqrt(sum(sum(gcint.I_std.^2, 3), 2) + (sum(s_gc.diode)/sum(s_air.diode)).^2 * sum(sum(airint.I_std.^2, 3), 2));

gc = load(gc_file);
qmin = max(q_gc(1), gc(1,1));
qmax = min(q_gc(end), gc(end,1));
qind = find(qmin < gc(:,1) & gc(:,1) < qmax);
q = gc(qind, 1);
tth = 2*asin(lambda * q / (4*pi));

% Cross section per q-bin. Factors for thickness, transmission and solid angle
xsection_scale = (1 / 10) * 1/gc_transmission * pixel_size_mm^2 / (4*pi*det_dist_mm^2);
gc_xsection = xsection_scale * gc(qind,2);
gc_xsection_err = xsection_scale * gc(qind,3);

% interpolate and scale w. angle dependent pixel solid angle and tilt
gc_exp = interp1(q_gc, I_gc, q) ./ (cos(tth).^3);
gc_exp_err = interp1(q_gc, Ierr_gc, q) ./ (cos(tth).^3);

s2 = gc_xsection_err.^2 + gc_exp_err.^2;
% Solve for scaling factor, weigh with combined variance^-1
wscale = sum(gc_xsection.*gc_exp./s2) / sum(gc_exp.^2 ./ s2);
fitchi = sum((gc_xsection - wscale*gc_exp).^2./(gc_xsection_err.^2 + (wscale*gc_exp_err).^2));
clf()
hold on
errband(q, gc_xsection, gc_xsection_err, 'r');
errband(q, wscale*gc_exp, wscale*gc_exp_err, 'b');
legend('Cross section', '    1 std', 'Experimental', '    1 std');
hold off

% scaling without error weighing
% scale = gc_exp \ gc_xsection;
%semilogy(q, gc_xsection, '*', q, scale*gc_exp)

% incoming flux (photons/s) determined for each q-channel:
%plot((1/gc_time) * gc_exp ./ gc_xsection)
inphotons = 1/gc_time * mean(gc_exp./gc_xsection);
diodescale = inphotons / (sum(s_gc.diode)/gc_time/gc_transmission);
fprintf('Glassy carbon transmission:        %g\n', gc_transmission);
fprintf('Chi^2 to known cross-section:      %g\n', fitchi);
fprintf('Flux on sample:                    %g ph/sec\n', inphotons);
fprintf('Scaling factor for diode readings: %g\n', diodescale);


function errband(x, y, yerr, colour);
%function errband(x, y, yerr, colour);
%
% Plot an error band between y-yerr and y+yerr.
% The colour defaults to blue.

if nargin<4
    colour='b';
end

x = x(:).';
y = y(:).';
yerr = yerr(:).';

lower = y-yerr;
upper = y+yerr;

hold on
plot(x, y, colour);
h = fill([x, fliplr(x)], [upper, fliplr(lower)], colour);
alpha(h, 0.5);
