% GET_ATT_LENGTH returns attenuation length of a material for a given energy (range)
%   formula...              chemical formula
%   energy...               single value in keV or energy range in keV
%   (optional) dens...      density, negative number for default values
%   (optional) ang...       grazing angle (default 90)
%   (optional) npts...      number of points
%   (optional) plot...      set to 1 for plotting
%
%   returns 
%       att...            (energy in keV, transmission)
%       req_density...      density in g/cm^3
%
%   examples:
%   get_att_length('Au', 8.7, -1, 45)
%   get_att_length('Pb', [11.2 24], 0.1, -1, 100)
%
%   03/2017

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

function [ att, req_density ] = get_att_length( formula, energy, varargin )

% check density
if nargin < 3
    dens = -1;
else
    dens = varargin{1};
end

% check angle
if nargin < 4
    ang = 90;
else
    ang = varargin{2};
end

% check specified number of points
if nargin < 5
    npts = 99;
else
    npts = varargin{3}-1;
end

% check if plotting is requested
if nargin < 6
    plot_att = false;
else
    plot_att = varargin{4};
end


% convert to keV
energy = energy * 1000;

if size(energy) ==1
    emin = energy-1;
    emax = energy+1;
    npts = 2;
    req_range = false;
elseif size(energy,2) == 2
    emin = energy(1);
    emax = energy(2);
    req_range = true;
else
    error('Only one specific energy or an energy range is supported.')
    
end

% check the energy range
if emin < 30 || emax > 30000
    error('Energies must be in the range 0.03 keV to 30 keV.')
end

% request the data
server = 'http://henke.lbl.gov/';
req = sprintf('Material=Enter+Formula&Formula=%s&Density=%f&Scan=Energy&Min=%d&Max=%d&Npts=%d&Fixed=%f&Plot=Log&Output=Plot', formula, dens, emin, emax, npts, ang);
data_req = webwrite('http://henke.lbl.gov/cgi-bin/atten.pl', req);

% find and read dat file
f_pos = strfind(data_req, '/tmp');
data_req = strsplit(data_req(f_pos(1):end), '.');
data = webread([server data_req{1} '.dat']);

% split data by line breaks
data = strsplit(data, '\n');

% extract density
req_density = data{1};
req_density = strsplit(req_density, '=');
req_density = strsplit(req_density{2}, ' ');
req_density = str2double(req_density{1});


% output
if npts==2 && ~req_range
    att = zeros(1,2);
    req_att = data{4};
    req_att = strsplit(req_att,' ');
    att(1,1) = str2double(req_att{2})/1000;
    att(1,2) = str2double(req_att{3});
else
    att = zeros(npts+1,2);
    for i=1:npts+1
        req_att = data{i+2};
        req_att = strsplit(req_att,' ');
        att(i,1) = str2double(req_att{2})/1000;
        att(i,2) = str2double(req_att{3});
    end
end

% plot attenuation
if plot_att
    if ~req_range
        fprintf('Requested plot for a single point.')
    end
    figure(76);
    plot(att(:,1), att(:,2))
    ylabel('attenuation length')
    xlabel('energy in keV')
    grid on;
    
end

end

