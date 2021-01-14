%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: 
%
% dose_calc(ptycho_recon, ptycho_data, param)
%
% Description:
%
% The function (1) takes one reconstruction and its data and estimate the
% dose (2) saves the dose estimation into a .txt file.
%
% Input: 
%
% ptycho_recon: reconstruction, including object, probe, and p
% ptycho_data: data for the reconstruction
% param_dose.mu = 1/(451*1e-6);     % 1/attenuation_length   in 1/m   (for CH2 @6.2keV)
%                                   % 1/(152.7*1e-6) for zeolite Na2Al2Si3O102H4O with 2 g/cm3 density at 6.2 keV
% param_dose.rho = 1000;            % Density in kg/m^3
% param_dose.setup_transmission = 0.55; % Intensity transmission of sample 
%                                       % (e.g. air path after the sample, windows, He, detector efficiency)
%                                       % 0.943 for 700 cm He gas at 760 Torr and 295 K @ 6.2 keV
%                                       % 0.780 for 10 cm air at 760 Torr and 295 K @ 6.2 keV
%                                       % 0.976 for 13 micron Kapton (polymide) with 1.43
%                                       % g/cm3 @ 6.2 keV
%                                       % 0.841 for 7 micron muskovite mica
%                                       % (KAl3Si3O11.8H1.8F0.2) with 2.76 g/cm3 @ 6.2 keV
%                                       % 0.914 for 5 cm of air at 6.2 keV 750 Torr 295 K
%                                       % 0.55 for 300 micron of mylar C10H8O4 with density 1.38 g/cm3 at 6.2 keV
% param_dose.overhead = 0.0;    % Extra dose during movement overhead, only applicable 
%                               % if shutter is not closed between exposures
% param_dose.fmask
% param_dose.scan_number 
% param_dose.num_proj 
% param_dose.output_folder (default: ./)
%
% Output:
%
% one jpg for one_data_frame
% one jpg for photons_per_shot_all
% one jpg for photons_per_obj_pix
% one txt for dose_estimate
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

function dose_calc(ptycho_recon, ptycho_data, param)
import utils.* 


%Estimating detected photons
probe = ptycho_recon.probe;
object = ptycho_recon.object;
p = ptycho_recon.p;
data = ptycho_data.data;
fmask = ptycho_data.fmask; 

if isfield(param,'mu') && isfield(param,'rho') && isfield(param,'setup_transmission') && isfield(param,'overhead') 
    mu = param.mu;
    rho = param.rho;
    setup_transmission = param.setup_transmission;
    overhead = param.overhead;
else
    error('Please specify param.mu, param.rho, param.setup_transmission, and param.overhead.\n');
end 

if isfield(param,'num_proj')
    num_proj = param.num_proj;
else
    verbose(0,'Using num_proj = 1');
    num_proj = 1;
end 

if isfield(param,'scan_number')
    scan_number = param.scan_number;
else
    scan_number = [];
end 

if isfield(param,'output_folder')
    output_folder = param.output_folder;
else
    verbose(0,'Using output_folder = ./');
    output_folder = '.';
end 

data = data .* fmask; 
photons_per_shot_all = sum(sum(data)); 
photons_per_shot = max(photons_per_shot_all);  

%Normalizing probe to photons per shot
probe_norm = sum(abs(probe).^2,3); 

probe_norm = probe_norm/sum(probe_norm(:));
probe_norm = probe_norm*photons_per_shot;
%
asize = size(probe);
%objaux = object*0;
illum_sum = zeros(size(object,1)+10,size(object,2)+10);


scanfirstindex = [1 cumsum(p.numpts)+1]; % First index for scan number
for ii = 1:p.numscans
    p.scanindexrange(ii,:) = [scanfirstindex(ii) scanfirstindex(ii+1)-1];
    p.scanidxs{ii} = p.scanindexrange(ii,1):p.scanindexrange(ii,end); 
end

for ii = p.scanidxs{1}
    Indy = round(p.positions(ii,1)) + [1:asize(1)];
    Indx = round(p.positions(ii,2)) + [1:asize(2)];
    illum_sum(Indy,Indx) = illum_sum(Indy,Indx)+probe_norm;
end

illum_sum = illum_sum(asize(1)/2:end-asize(1)/2,asize(2)/2:end-asize(2)/2);

% in case of laminography the field of view isnot rectangular ->
% exclude the empty regions in the illumination function 
illum_mask = illum_sum > mean(illum_sum) * 0.1; 

flux_in_area = sum(sum(illum_sum .* illum_mask)); %photons
area = sum(illum_mask(:))*p.dx_spec(1)^2;  % meters^2

I = flux_in_area/area;  %ph/meters^2
hv = 9.9334947e-16*(p.energy/6.2);      %6.2keV in joules

D = mu*I*hv*num_proj/rho;
D_with_gas = D/setup_transmission;
D_with_overhead = D_with_gas*(1+overhead);
verbose(0,'**********************************************')
verbose(0,'Dose report for %d projections, Scan %d',num_proj, scan_number)
verbose(0,'**********************************************')
verbose(0,'Measured photons per frame = %.2e photons',photons_per_shot);
verbose(0,'N_0 used for imaging for one projection = %.2e photons/micron^2',I*1e-12);
verbose(0,'Dose used for imaging, D = %.2e Gy',D)
verbose(0,'Accounting for experiment transmission, D = %.2e Gy',D_with_gas)
verbose(0,'Accounting for experiment transmission and overhead, D = %.2e Gy',D_with_overhead)

%=======================
figure(1); clf
plotting.imagesc3D(log10(1+data));
caxis([0,log10(max(data(:)))])
colormap(plotting.franzmap); colorbar
axis xy equal tight
title('Data frames, log10')
filename = fullfile(output_folder,sprintf('/%s_one_data_frame.jpg',p.run_name));
verbose(1,'saving %s',filename);
print('-djpeg','-r300',filename);

figure(2); clf
plot(squeeze(photons_per_shot_all)); grid on;
title(['Number of measured photons per frame = ' num2str(photons_per_shot)]); 
filename = fullfile(output_folder,sprintf('/%s_photons_per_shot_all.jpg',p.run_name));
verbose(1,'saving %s',filename);
print('-djpeg','-r300',filename);

figure(3); clf
imagesc(illum_sum)
colormap(plotting.franzmap); colorbar
axis image xy 
title('Photons per pixel of the object')
filename = fullfile(output_folder,sprintf('/%s_photons_per_obj_pix.jpg',p.run_name));
verbose(1,'saving %s',filename);
print('-djpeg','-r300',filename);

%=======================
filename = fullfile(output_folder,sprintf('/%s_dose_estimate_S%05d.txt',p.run_name, scan_number));
fid = fopen(filename,'w');
fprintf(fid,'Scan = %d\n',scan_number);
fprintf(fid,'Measured photons per frame = %.3e\n',photons_per_shot);
fprintf(fid,'N_0 used for imaging for one projection (I*1e-12) = %.3e photons/micron^2\n',I*1e-12);
fprintf(fid,'num_proj = %d\n',num_proj);
fprintf(fid,'hv = %.5e\n\n',hv);
fprintf(fid,'D = mu*I*hv*num_proj/rho \n');
fprintf(fid,'D_with_gas = D/setup_transmission \n');
fprintf(fid,'D_with_overhead = D_with_gas*(1+overhead) \n\n');
fprintf(fid,'If mu = %.3e m^-1, rho = %.3e kg/m^3, setup_transmission = %.3f, then:\n', mu, rho, setup_transmission);
fprintf(fid,'Accounting for experiment transmission, D_with_gas = %.3e Gy = %.3f MGy\n\n', D_with_gas, D_with_gas/1e6);
fprintf(fid,'asize = %d, pixel size = %.4f nm\n', asize(1), p.dx_spec(1)*1e9);
fclose(fid);