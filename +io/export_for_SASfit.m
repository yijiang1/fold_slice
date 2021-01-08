% This script is to plot, correct and export solution SAXS data to SASfit
% accounts for transmission, time and thickness correction
% scales the data to a calibration factor 
% background correction, removal of bad pixels
% not suitable for anisotropic data
% saves the output to be used in SASfit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  EDIT HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% give the calibration factor for absolute intensity, calculated previously
cal_factor_SAXS = 2.91e-4;
cal_factor_WAXS = 2.01e-5 ;

% where the data is saved (Data10, afs, p-account)
base_dir = '~/Data10/';%'/sls/X12SA/Data20/e16598/';
save_dir = '~/Data10/';%'/mnt/das-gpfs/work/p16598/';
eaccount = beamline.identify_eaccount; % 'e16598';

% samples and thicknesses
Air = 15; % scan used for transmission calculation
sample_scan = [34:38]; % should be given
sample_thickness = 0.15; % in cm: important for absolute scattering
back_scan = []; % used as background, leave it empty [] for no subtraction !!NOT TESTED!!
back_thickness = 0.01; % in cm: important for absolute scattering

% export data for SASfit?
export_sasfit = 1;

%plot curves?
plot_curves = 0;

% save figures?
save_fig = 0;

% average the scan points? 1 = yes, 0 = no
average_scan = 1;
% scale also the waxs data? yes = 1; no = 0;
use_waxs = 0;

% which bad pixels should be removed
bad_pixel = []; %given as a vector [811, 825]

% use all data measurement points or skip some (faster)
skip_measurements = [100]; % use 1 to show all

% used to reduce noise at the beginning and end of scattering curve
skip_first_points = 55;
skip_last_points = 150;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load the diode value for the air
S_air = io.spec_read(base_dir,'ScanNr',Air);
%scale in case the exposure times are different
exp_time = S_air.sec(1,1);
scale_air = 1/exp_time;
Air_data = load(sprintf('%s/analysis/radial_integration/%s_1_%05d_00000_00000_integ.mat', base_dir, eaccount, Air));
I_air = mean(Air_data.I_all, 3);
% average over the segments when needed
if size(I_air, 2) > 1
    I_air = (I_air .* Air_data.norm_sum)./sum(Air_data.norm_sum, 2);
    I_air = sum(I_air, 2);
end

if use_waxs
    Air_waxs = load(sprintf('%s/analysis/radial_integration_waxs/%s_2_%05d_00000_00000_integ.mat', base_dir, eaccount, Air));
    I_air_waxs = mean(squeeze(Air_waxs.I_all), 2);
end
%% background correction
if ~isempty(back_scan)
    for b = 1:length(back_scan)
        bgr = load(sprintf('%s/analysis/radial_integration/%s_1_%05d_00000_00000_integ.mat', base_dir, eaccount, back_scan(b)));
        S_back = spec_read(base_dir,'ScanNr',back_scan(b));
        % in case the burst scan takes place, the transmission is
        % calculated differently
        if ~isempty(findstr(S_back.S, 'burst_scan'))
            delimiter = ' ';
            formatSpec = '%*s%*s%s%[^\n\r]';
            fileID = fopen(sprintf('%smcs/S00000-00999/S%05d/%s_%05d.dat', base_dir,back_scan(b), eaccount, back_scan(b)), 'r');
            dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);
            exp_time = dataArray{1,1}{7,1};
            scale_back = 1/str2num(exp_time);
            diode = mean(str2num(dataArray{1,1}{9,end}));
            transm_back = ((diode*scale_back)/(mean(S_air.diode)*scale_air));
            fclose(fileID);
        else
            %scale in case the exposure times are different
            exp_time = S_back.sec(1,1);
            scale_back = 1/exp_time;
            transm_back = (mean(S_back.diode)/mean(S_back.bpm4i))/(mean(S_air.diode)/mean(S_air.bpm4i));
        end
        %average background
        I_bgr = mean(bgr.I_all, 3);
        if size(I_bgr, 2) > 1
            I_bgr = (I_bgr .* bgr.norm_sum)./sum(bgr.norm_sum, 2);
            I_bgr = sum(I_bgr, 2);
        end
        I_bgr = ((((I_bgr*scale_back)*1/transm_back)-(I_air*scale_air))*1/back_thickness);
        I_bgr = I_bgr * cal_factor_SAXS;
        if use_waxs
            %average background_WAXS
            bgr_waxs = importdata(sprintf('%s/analysis/radial_integration_waxs/%s_2_%05d_00000_00000_integ.mat', base_dir, eaccount, back_scan(b)));
            I_bgr_waxs = mean(squeeze(bgr_waxs.I_all), 2);
            I_bgr_waxs = ((((I_bgr_waxs*scale_back)*1/transm_back)-(I_air_waxs*scale_air))*1/back_thickness);
            I_bgr_waxs = I_bgr_waxs * cal_factor_WAXS;
        end
    end
else
    I_bgr_waxs = 0;
    I_bgr = 0;
end

%% load and correct the sample
for s = 1:length(sample_scan)
    sample_filename=sprintf('%s/analysis/radial_integration/%s_1_%05d_00000_00000_integ.mat', base_dir, eaccount, sample_scan(s));
    if exist(sample_filename) == 2
        display(['reading file ',sample_filename])
        sample = load(sample_filename);
    else
        continue
    end
    S_s = io.spec_read(base_dir,'ScanNr',sample_scan(s));
    if ~isempty(findstr(S_s.S, 'burst_scan'))
        delimiter = ' ';
        formatSpec = '%*s%*s%s%[^\n\r]';
        fileID = fopen(sprintf('%smcs/S00000-00999/S%05d/%s_%05d.dat', base_dir,sample_scan(s), eaccount,sample_scan(s)), 'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);
        exp_time = dataArray{1,1}{7,1};
        scale_s = 1/str2num(exp_time);
        diode = mean(str2num(dataArray{1,1}{9,end}));
        transm_sample = ((diode*scale_s)/(mean(S_air.diode)*scale_air));
        fclose(fileID);
    else
        %scale in case the exposure times are different
        exp_time = S_s.sec(1,1);
        scale_s = 1/exp_time;
        transm_sample = (mean(S_s.diode)/mean(S_s.bpm4i))/(mean(S_air.diode)/mean(S_air.bpm4i));
    end
    %load the sample
    I_sample = squeeze(sample.I_all);
    q_sample = sample.q';
    
    if use_waxs
        %average background_WAXS
        sample_waxs = load(sprintf('%s/analysis/radial_integration_waxs/%s_2_%05d_00000_00000_integ.mat', base_dir, eaccount, sample_scan(s)));
        I_sample_waxs = (sample_waxs.I_all);
        q_sample_waxs = sample_waxs.q';
    else
        q_sample_waxs = [];
        I_sample_waxs = [];
    end
    if average_scan
        I_sample = mean(I_sample, 3);
        I_std = mean(sample.I_std, 3);
        if size(I_sample, 2) > 1
            I_sample = (I_sample .* sample.norm_sum)./sum(sample.norm_sum, 2);
            I_std = (I_std .* sample.norm_sum)./sum(sample.norm_sum, 2);
            I_std = sum(I_std, 2).*cal_factor_SAXS;
            I_sample = sum(I_sample, 2);
        end
        I_std = I_std(skip_first_points:end-skip_last_points,:);
        I_sample = ((((I_sample*scale_s)*1/transm_sample)-(I_air*scale_air))*1/sample_thickness);
        if ~isempty(bad_pixel)
            I_sample(bad_pixel,1) = (I_sample(bad_pixel-1,1)+I_sample(bad_pixel+1,1))/2;
        end
        I_sample = I_sample * cal_factor_SAXS;
        I_cor = (I_sample-I_bgr);
        I_cor = I_cor(skip_first_points:end-skip_last_points,:);
        
        if use_waxs
            hold on
            I_sample_waxs = median(I_sample_waxs,3);
            I_sample_waxs = ((((I_sample_waxs.*scale_s).*1/transm_sample)-(I_air_waxs.*scale_air)).*1/sample_thickness);
            I_sample_waxs = I_sample_waxs * cal_factor_WAXS;
            I_cor_waxs = (I_sample_waxs-I_bgr_waxs);
            
        else
            I_cor_waxs = [];
        end
        I_total = [I_cor; I_cor_waxs];
        
        q_total = [q_sample(skip_first_points: end-skip_last_points,:); q_sample_waxs];
        [q_total, index] = sort(q_total);
        I_total = I_total(index);
        
        if plot_curves
            figure
            plot(q_total*10, I_total);
            set(gca,'XScale','log', 'YScale','log');
            grid on;
            box on;
            xlabel('scattering vector q (nm^{-1})');
            ylabel('differential scattering cross-section (cm^{-1})');
            hold on
        end
        if export_sasfit
            save_data = [q_total*10, I_total, I_std];
            filename = sprintf('scan_%05d_avg', sample_scan);
            save(sprintf('%sanalysis/dat_files/%s.dat', save_dir , filename) , 'save_data', '-ascii');
        end
    else
        if plot_curves
            figure
            hold on
        end
        for i = 1:skip_measurements:size(sample.I_all, 3)
            I_point = sample.I_all(:,:,i);
            I_point_std = sample.I_std(:,:, i);
            if size(I_point, 2) > 1
                I_point = (I_point .* sample.norm_sum)./sum(sample.norm_sum, 2);
                I_point = sum(I_point, 2);
                I_point_std = (I_point_std .* sample.norm_sum)./sum(sample.norm_sum, 2);
                I_point_std = sum(I_point_std, 2).*cal_factor_SAXS;
            end
            I_point_std = I_point_std(skip_first_points:end-skip_last_points,:);
            if ~isempty(bad_pixel)
                I_point(bad_pixel,1) = (I_point(bad_pixel-1,1) + I_point(bad_pixel+1,1))/2;
            end
            I_point = ((((I_point*scale_s)*1/transm_sample)-(I_air*scale_air))*1/sample_thickness);
            I_point = I_point * cal_factor_SAXS;
            I_cor = (I_point-I_bgr);
            I_cor = I_cor(skip_first_points: end-skip_last_points,:);
            
            if use_waxs
                I_point_waxs = I_sample_waxs(:,i);
                I_point_waxs = ((((I_point_waxs*scale_s)*1/transm_sample)-(I_air_waxs*scale_air))*1/sample_thickness);
                I_point_waxs = I_point_waxs * cal_factor_WAXS;
                I_cor_waxs = (I_point_waxs-I_bgr_waxs);
                I_point_std_WAXS = I_sample_waxs(:,:, i);
                I_point_std_WAXS = I_point_std_WAXS.*cal_factor_WAXS;
                I_point_std_WAXS = sum(I_point_std_WAXS, 2).*cal_factor_SAXS;
            else
                I_cor_waxs = [];
            end
            I_total = [I_cor; I_cor_waxs];
            q_total = [q_sample(skip_first_points: end-skip_last_points,:); q_sample_waxs];
            [q_total, index] = sort(q_total);
            I_total = I_total(index);
            I_point_std_total= [I_point_std; I_point_std_WAXS];
            if plot_curves
                plot(q_total*10, I_total);
                grid on;
                box on;
                set(gca,'XScale','log', 'YScale','log');
                xlabel('scattering vector q (nm^{-1})');
                ylabel('differential scattering cross-section (cm^{-1})');
                axis tight
                hold on
                drawnow
            end
            if export_sasfit
                save_data = [q_total*10, I_total, I_point_std_total];
                filename = sprintf('scan_%05d_pt_%05d', sample_scan(s), i);
                save(sprintf('%sanalysis/dat-files/%s.dat', save_dir , filename) , 'save_data', '-ascii');
            end
        end
        
    end
    
end
if save_fig
    %save the results
    saveas(gcf, sprintf('%sanalysis/scanNr_%05d.jpg', save_dir , sample_scan))
end

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
