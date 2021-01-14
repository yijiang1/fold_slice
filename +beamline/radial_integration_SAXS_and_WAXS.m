% radial_integration_SAXS_and_WAXS.m
% Template for radial integration made around 2015
% Changes:
% 2016-08-22: define mask files at the beginning, allowing for a flag in case it needs to be repeated  
% add the save fast and v6

% License at the end of script

clear all
close all

%% step 0: add the path for the matlab-scripts (fill in userID,detno and specdatfile)
addpath ..
%e-account followed by underline
userID   = [beamline.identify_eaccount '_']; 
% detector number: 1 for SAXS Pilatus 2M, 2 for WAXS Pilatus 300k, 3 for
% SAXS Eiger 500 k
detno = 1;
% which data format to save? '-v6' is the standard.
save_format = '-v6';
% flag for filenames for valid pixel mask, beamstop mask coordinates and integration mask. 
% Example: '_2M_at_two_meters'
% Leave empty '' for default folder and filenames. 
file_flag='';  

% change here for offline analysis
homedir  = sprintf('~/Data10/');
%homedir  = '/mnt/das-gpfs/work/p16268/';

%CHANGE: spec dat file
SpecDatFile = '~/Data10'; 

if (detno == 1 )||(detno == 2)
    datadir  = fullfile(sprintf('%s',homedir),sprintf('pilatus_%d/',detno));
elseif detno == 3
    datadir  = fullfile(sprintf('%s',homedir),sprintf('eiger'));
end
if detno == 2
    integdir = sprintf('%sanalysis/radial_integration_waxs%s/',homedir,file_flag);
elseif detno == 1
    integdir = sprintf('%sanalysis/radial_integration%s/',homedir,file_flag);
elseif detno == 3
    integdir = sprintf('%sanalysis/radial_integration_eiger%s/',homedir,file_flag);
end
if detno == 2
    outdir   = sprintf('%sanalysis/data_waxs%s/',homedir,file_flag);
elseif detno == 1
    outdir   = sprintf('%sanalysis/data/%s',homedir,file_flag);
elseif detno == 3
    outdir   = sprintf('%sanalysis/data_eiger%s/',homedir,file_flag);
end
addpath(sprintf('%smatlab/',homedir));
if (detno == 1 )||(detno == 2)
    maskfilename = sprintf('%spilatus_%d_valid_mask%s.mat', outdir,detno,file_flag); 
    integmaskfilename=sprintf('%spilatus_%d_integration_masks%s.mat',outdir,detno,file_flag); 
elseif detno == 3
    maskfilename = sprintf('%seiger_%d_valid_mask%s.mat', outdir,detno,file_flag); 
    integmaskfilename=sprintf('%seiger_%d_integration_masks%s.mat',outdir,detno,file_flag); 
end
maskcoordfilename=sprintf('%smask_coordinates_%d%s.mat',outdir,detno, file_flag); 

dirs = whos('-regexp','.*dir$');
for ii=1:numel(dirs)
  dir_to_do = eval(dirs(ii).name);
  if ~exist(dir_to_do,'dir')
    fprintf('creating directory %s\n', dir_to_do);
    system(sprintf('mkdir -p %s',dir_to_do));
  end
end
%% enter scan numbers of standards
%glassy carbon, glassy carbon moved detector to side, air scattering, first
%one is glassy carbon used to remove beamstop later
 scannr = [14 15 14];
 %AgBE (for SAXS and WAXS), LaB6 (for WAXS), Si (for WAXS)
   todo      = [12 13 14];
  legendstr = {'AgBE';'LaB6';'Si'};

%% step 1: prepare the valid pixel mask
redo = 1;

if (redo)
  fprintf('preparing the valid pixel mask\n');

% calculating the union of several valid pixel masks
%   starting with a rather dark file to discriminate hot pixels
  system(sprintf('rm -f %s', maskfilename));

if (detno == 1 )||(detno == 2)
    prepvalidmask_args   = {};
    compilex12sa_args    = {'DetectorNumber',detno,'FileExtension','cbf'};
    integrate_range_args = {'PilatusDetNo',detno,'FileExtension','cbf'};
elseif detno == 3
    prepvalidmask_args = {'H5Location','/eh5/images/','FilenameMask','*'};
    compilex12sa_args = {'FileExtension','h5'};
end
  
  for ii=scannr
    beamline.prep_valid_mask(utils.compile_x12sa_filename(ii,-1, ...
        'BasePath',datadir,'BaseName',userID,compilex12sa_args{:}), ...
        'ThresholdDark',1, ...
        'ThresholdHot',20, ...
        'Extend','or', ...
        'FilenameValidMask',maskfilename,prepvalidmask_args{:});
      %  'FigNo',ii==scannr(end));
  end  
end
%% step 2: cut out beam stop and shadows manually (for WAXS only necessary if there is a shadow)

redo = 1;
if (redo)
  scannr = scannr(1);
  if (detno == 1)||(detno == 2)
      compilex12sa_args = {'DetectorNumber',detno,'FileExtension','cbf'};
      imageshow_args = {};
  elseif (detno == 3)
      compilex12sa_args = {'FileExtension','h5'};
      imageshow_args = {'H5Location','/eh5/images/'};
  end
  % include the beamstop in the valid pixel mask - follow instructions in
  % popup box
  beamline.choose_beamstop_mask(utils.compile_x12sa_filename(scannr(1),0, 'BasePath',datadir,'BaseName',userID, compilex12sa_args{:}),...
      'ReadCoord',0,'SaveCoord',1, 'SaveData',1,'FilenameValidMask',maskfilename,'FilenameCoord',maskcoordfilename, 'ImageShowArgs', imageshow_args)
  
end
%% show silver behenate scattering to find the radius of the first ring (only SAXS)
if (detno==1)
plotting.image_show(utils.compile_x12sa_filename(todo(1),0, ...
          'PointWildcard', 1, ...
          'SubExpWildcard', 1, ...
          'DetectorNumber',detno, ...
          'BasePath',datadir,'BaseName',userID), ...
          'IsFmask', true);
elseif (detno == 3)
    filepath = utils.compile_x12sa_dirname(todo(1));
    D = dir(fullfile(datadir,filepath,'*.h5'));
    plotting.image_show(fullfile(D(1).folder,D(1).name), ...
          'H5Location','/eh5/images/');
end
%% here you have to give some manual inputs to run step 3
% for SAXS you have to put y pixel value of the the silver behenate ring above the beamstop, and the order of the peak that you chose
if (detno==1)||(detno == 3)
    order_AgBE = 1;
    y_from = 509;
    y_to   = 514;
    cen_guess = [];  %[y,x] ; leave empty, i.e. cen_guess=[], for automatic guess;
    %and choose how many sectors you want to do the integration (16 for
    %anisotropic scattering, 1 for isotropic scattering
    num_segments=16;
elseif (detno==2)
    %for WAXS you can run with the default values to start with and adjust in
    %case an error appears or the fit (shown in figure 4) is bad
    
    open('+beamline/WAXS_standards.fig');
    %give the order of the first silver behenate ring appearing
    %(compare with WAXS_standards.fig)
    order_AgBe=7;
    
    %parameter used in finding the x-position, default 5, if in figure 20 the
    %blue curve is all zeros, lower this value (necessary for low intensity of
    %silver behenate measurement
    d = 5;
    
    %threshold to find WAXS peak of standards, default is 50, might be lowered
    %for lower intensities
    threshold=[2 50 100];
    %if wrong peaks are found tune finding the right peaks with the window
    %where peaks are being searched here, default is min=0 and max=1500,
    %(see WAXS_standards.fig)
    min_AgBE=0;
    max_AgBE=1500;
    min_Si=0;
    max_Si=1500;
    min_LaB6=0;
    max_LaB6=1500;
    
end
% step 3: prepare integration mask
% For the WAXS mask this is still a bit clunky. You can adjust above the
% min and max values where it will look for a peak and the threshold. Also
% in the fit for the horizonal position make sure there is both red and
% blue peaks for the fitting, if not you can adjust the d parameter above.
% Decreasing it helps when the silver behenate scattering is low.

if (detno==1)
    scannr = todo(1);
else
    %here enter the scannumbers of the standards
    %   todo      = [211,208,212];
    %   legendstr = {'AgBE';'LaB6';'Si'};
    scannr    = todo(1);
    S = io.spec_read(SpecDatFile,'ScanNr',todo(1));
end

if (detno==1)||(detno==2)
    I = plotting.image_show(utils.compile_x12sa_filename(scannr,0, ...
          'PointWildcard', 1, ...
          'SubExpWildcard', 1, ...
          'DetectorNumber',detno, ...
          'BasePath',datadir,'BaseName',userID), ...
          'IsFmask', true);
elseif (detno == 3)
    filepath = utils.compile_x12sa_dirname(scannr);
    D = dir(fullfile(datadir,filepath,'*.h5'));
    I = plotting.image_show(fullfile(D(1).folder,D(1).name), ...
          'H5Location','/eh5/images/');
end


mask = getfield(load(maskfilename),'valid_mask');
mask.frame = zeros(mask.framesize);
mask.frame(mask.indices) = 1;

I = mean(I.data,3).*mask.frame;
if (detno==1)||(detno == 3)
    J = ifftn(fftn(I,size(I)*2-[1 1]).^2);
    if isempty(cen_guess)
        cen_guess = math.peakfit2d(J)/2; %[y,x]
    end
    
    if (detno == 1)
        filename_center = utils.compile_x12sa_filename(scannr(1),0, 'BasePath',datadir,'BaseName',userID);
        imageshow_args = {};
    elseif (detno == 3)
        filename_center =  fullfile(D(1).folder,D(1).name);
        imageshow_args =  {'H5Location','/eh5/images/'};
    end
    
    [cen]=utils.get_beam_center(filename_center,'GuessX',cen_guess(2),'GuessY',cen_guess(1), ...
        'RadiusFrom',y_from-cen_guess(1),'RadiusTo',y_to-cen_guess(1), ...
        'TestX',4,'TestY',4,'FilenameValidMask',maskfilename, imageshow_args{:});
    
else
    % this isn't nice yet
    % i)   it depends on the chosen orientation on how to read
    %      detector-2 images
    % ii)  it merely finds maximum values instead of fitting, possibly
    %      with sub-pixel precision
    % iii) as a consequence, figuring out which values are trustworthy
    %      is done rather crudly
    %d = 3; %5 seams not to work if intensity of silver behenate is too low??
    if (detno == 2)
        imageshow_args = {};
    end
    dx = 30;
    [s1,s2] = size(I);
    
    J = ifft(fft(I,s1*2-1,1).^2,[],1);
    [~,n] = max(J);
    w = std(I,1,1)./sqrt(mean(I,1));
    o = 1:numel(n);
    
    o = o(w>d);
    n = n(w>d)/2;
    
    o = o(abs(n-s1/2)<dx);
    n = n(abs(n-s1/2)<dx);
    
    x = s1/2+linspace(-dx,dx,4*dx+1);
    
    figure(20)
    m = histc(n,x);
    [~,n0] = max(m);
    plot(x,m)
    hold on
    
    s = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[  0,s1/2-dx,  0,  0,  0],...
        'Upper',[Inf,s1/2+dx,Inf,Inf,Inf],...
        'Startpoint',[10,x(n0),1,10,1]);
    f = fittype('a*exp(-((x-b)/c)^2)+d*exp(-((x-n)/e)^2)', ...
        'problem','n','options',s);
    [c,~] = fit(x',m',f,'problem',s1/2);
    figure(50)
    plot(c,'r');
    hold off
    
    figure(10)
    cen1 = c.b;
    o = o(abs(n-cen1)<=1);
    n = n(abs(n-cen1)<=1);
    hold on
    plot(o,n,'w.')
    plot([1 s2],[1 1]*round(cen1),'w')
    x = 1:s2;
    plot(x(mask.frame(round(cen1),:)>0), ...
        log(I(round(cen1),mask.frame(round(cen1),:)>0))/ ...
        max(log(I(round(cen1),mask.frame(round(cen1),:)>0)))*s1, ...
        'k')
    hold off
    
    figure(30)
    WAXS = zeros(s2,numel(todo));
    WAXS(:,1) = I(round(cen1),:);
    for ii=2:numel(todo)
        I = io.image_read(utils.compile_x12sa_filename(todo(ii),0, ...
            'PointWildcard', 1, ...
            'SubExpWildcard', 1, ...
            'DetectorNumber',detno, ...
            'BasePath',datadir,'BaseName',userID), ...
            'IsFmask', 1);
        WAXS(:,ii) = mean(I.data(round(cen1),:,:),3);
    end
    
    h = semilogy(WAXS);
    legend(legendstr)
    
    % finding peaks "automatically"
    x_coord = [];
    q_coord = [];
    hold on
    peaks = cell(1,size(WAXS,2));
    for ii=1:size(WAXS,2)
        %the treshhold value, default set to 50, might be adjusted
        peaks{ii} = utils.peakfinder((WAXS(:,ii)),threshold(ii));
        %peaks{ii} = peakfinder((WAXS(:,ii)),50);
        if strcmp(legendstr{ii},'AgBE')
            tmp = peaks{ii};
            tmp = tmp(tmp>=min_AgBE);
            peaks{ii} = tmp(tmp<=max_AgBE);
            
        end
        if strcmp(legendstr{ii},'Si')
            tmp = peaks{ii};
            tmp = tmp(tmp>=min_Si);
            peaks{ii} = tmp(tmp<=max_Si);
        end
        if strcmp(legendstr{ii},'LaB6')
            tmp = peaks{ii};
            tmp = tmp(tmp<=max_LaB6);
            peaks{ii} = tmp(tmp>=min_LaB6);
            
        end
        
        
        x_coord = vertcat(x_coord,peaks{ii});
        if strcmp(legendstr{ii},'AgBE')
            q0 = 2*pi/58.38;
            q_coord = horzcat(q_coord,q0*(order_AgBe+(0:numel(peaks{ii})-1)));
        elseif strcmp(legendstr{ii},'LaB6')
            q0 = 2*pi/4.1549;
            q_coord = horzcat(q_coord,q0*sqrt((1:numel(peaks{ii}))));
        elseif strcmp(legendstr{ii},'Si')
            q0 = 2*pi/5.4308;
            q_coord = horzcat(q_coord,q0*sqrt(3));
        end
        semilogy(peaks{ii},WAXS(peaks{ii},ii),'.', ...
            'Color',get(h(ii),'Color'), ...
            'MarkerSize',24)
    end
    hold off
    
    figure(40); clf
    if (numel(x_coord)>3)
        %       fprintf('%f\t%f\n',[x_coord';q_coord])
        %       % a
        %       % b
        %       % c
        s = fitoptions('Method','NonlinearLeastSquares',...
            'Lower'     ,[-Inf,-Inf,  0],...
            'Upper'     ,[ Inf,   0,1e3],...
            'Startpoint',[s2/2, 200,550]);
        f = fittype('4*pi/l*sin((atan((a-b)*p/c)+atan((x-a)*p/c))/2)', ...
            'problem',{'p','l'},'options',s);
        [c,~] = fit(x_coord,q_coord',f,'problem',{.172,12.398/S.mokev});
        subplot(2,1,1)
        plot(x_coord,q_coord,'x');
        hold on
        drawnow;
        tmp = axis;
        x = linspace(c.b,tmp(2));
        plot(x,feval(c,x),'r');
        subplot(2,1,2)
        bar(x_coord,feval(c,x_coord)-q_coord');
        xlim(tmp(1:2));
        dc = confint(c);
        dc = (dc(2,:)-dc(1,:))/2;
        fprintf(['detector distance:\t%.1fmm,    \t%.1fmm\n', ...
            'center of rings:  \t%.1fpixels,\t%.1fpixels\n', ...
            'angle of detector:\t%.1fdeg,   \t%.1fdeg.\n'],   ...
            c.c,dc(3), ...
            c.b,dc(2), ...
            atan((c.a-c.b)*c.p/c.c)/pi*180, ...
            180/pi*c.p/c.c*sqrt(dc(1)^2+dc(2)^2 + ((c.a-c.b)/c.c*dc(3))^2));
    end
end


tic
if (detno==1)||(detno==3)
    S = io.spec_read(SpecDatFile,'ScanNr',todo(1));
    fprintf('preparing the integration mask(s)\n');
    beamline.prep_integ_masks(utils.compile_x12sa_filename(todo(1),0, ...
        'BasePath',datadir,'BaseName',userID, compilex12sa_args{:}), ...
        cen, ...
        'DetNo',detno, ...
        'NoOfSegments',num_segments, ...
        'FilenameValidMask',maskfilename, ...
        'FilenameIntegMasks',integmaskfilename, imageshow_args{:});
    
    beamline.integrate_range(todo(1),todo(1),1, ...  % change for not re-running on already integrated files
        'OutdirData',integdir, ...
        'BasePath',datadir,'BaseName',userID, ...
        'FilenameIntegMasks',integmaskfilename, ...
        compilex12sa_args{:},imageshow_args{:});
    
elseif (detno==2)
    S = io.spec_read(SpecDatFile,'ScanNr',todo(1));
    fprintf('preparing the integration mask(s)\n');
    beamline.prep_integ_masks(utils.compile_x12sa_filename(todo(1),0, ...
        'DetectorNumber',detno, ...
        'BasePath',datadir,'BaseName',userID), ...
        [c.b cen1], ...
        'DetNo',detno, ...
        'Wavelength_nm', 12.398/S.mokev, ...
        'NormalXY', [c.a cen1], ...
        'DetDist_mm', c.c, ...
        'PixelSize_mm', .172, ...
        'NoOfSegments',1, ...
        'FilenameValidMask',maskfilename, ...
        'FilenameIntegMasks',integmaskfilename, ...
        'DisplayValidMask',0);
end
toc


%% calculate detector distance (SAXS only) check in Figure 100 if the peak_agbe really is the 1st order AgBE
if (detno==1)||(detno==3)
    [x,y] = plotting.plot_radial_integ(sprintf('%s%s%d_%05d_00000_00000_integ.mat',integdir,userID,1,todo(1)));
    %%the 1st order silver behenate is at ... pixels
    %peakfinder(log(y(10:end)),1);
    peaks2 = utils.peakfinder(log(y(10:end)),1);
    peak_agbe = x(peaks2(order_AgBE+1))+9 %normally the 1st order AgBE, check!
    wavelength = 12.398/S.mokev;
    detector_distance = peak_agbe*.172/tan(2*asin(wavelength*order_AgBE/(2*58.38)))
end
%% redo SAXS integration mask now it will take the detector distance into account and also save the q-value
if (detno==1)||(detno==3)
    if (detno == 1)
        detector_pixelsize = 0.172;
    elseif (detno == 3)
        detector_pixelsize = 0.075;
    end
    S = io.spec_read(SpecDatFile,'ScanNr',todo(1));
    fprintf('preparing the integration mask(s)\n');
    beamline.prep_integ_masks(utils.compile_x12sa_filename(todo(1),0, ...
        'BasePath',datadir,'BaseName',userID,compilex12sa_args{:}), ...
        cen, ...
        'DetNo',detno, ...
        'NoOfSegments',num_segments, ...
        'Wavelength_nm', 12.398/S.mokev, ...
        'DetDist_mm', detector_distance, ...
        'PixelSize_mm', detector_pixelsize, ...
        'FilenameValidMask',maskfilename, ...
        'FilenameIntegMasks',integmaskfilename, imageshow_args{:});
end
%% step 5: radial integration & averaging of files -- 
%start here again if you merely want to integreat
%for fast measurements (i.e. scanning SAXS) start on several cn parallel
%adjust therefor integrate_range(scan_no_from,scan_no_to,scan_no_step) 
%and rund only step 0 and step 5
save_format = '-v6';

close all
% beamline.integrate_range(107,1e8,3, ...  % change for not re-running on already integrated files
%     'PilatusDetNo',detno, ...
%     'OutdirData',integdir, ...
%     'BasePath',datadir,'BaseName',userID, ...
%     'FilenameIntegMasks',integmaskfilename, 'SaveFormat', save_format);

beamline.integrate_range(136,137,1, ...  % change for not re-running on already integrated files
    'OutdirData',integdir, ...
    'BasePath',datadir,'BaseName',userID, ...
    'FilenameIntegMasks',integmaskfilename, 'SaveFormat', save_format, ...
    integrate_range_args{:},imageshow_args{:});


%% or alternatively when computers node are ready and matlab is open
save_format = '-v6';
fprintf('beamline.integrate_range(107,1e8,4,''OutdirData'',''%s'',''BasePath'',''%s'',''BaseName'',''%s'',''FilenameIntegMasks'',''%s'',''SaveFormat'', ''%s''',integdir,datadir,userID,integmaskfilename,save_format)
args={'OutdirData', integdir,'BasePath',datadir ,'BaseName',userID ,'FilenameIntegMasks',integmaskfilename ,'SaveFormat',save_format };

for ii = 1:2:numel(integrate_range_args)
    if ischar(integrate_range_args{ii+1})
        straux = '''%s''';
    elseif isnumeric(integrate_range_args{ii+1})
        straux = '%d';
    end
    fprintf( [',''%s'',' straux ' '] ,integrate_range_args{ii},integrate_range_args{ii+1});
    args=[args,integrate_range_args{ii},integrate_range_args{ii+1}];
end
for ii = 1:2:numel(imageshow_args)
    if ischar(imageshow_args{ii+1})
        straux = '''%s''';
    elseif isnumeric(imageshow_args{ii+1})
        straux = '%d';
    end
    fprintf([',''%s'',' straux ' '],imageshow_args{ii},imageshow_args{ii+1});
    args=[args,imageshow_args{ii},imageshow_args{ii+1}];
end

% if detno==2 % Disable CReader for WAXS detector since currently it's not supported.
%     fprintf([',''CReader'',0 ']);
%     args=[args,'CReader',0];
% end

fprintf(');\n')

folder_todo=utils.abspath('~/Data10/analysis/radial_integration_todo/');
if ~exist(folder_todo)
    mkdir(folder_todo);
end

save(sprintf([folder_todo 'vargin_det%d.mat'],detno),'args');
fprintf(['Parameters saved to' folder_todo 'vargin_det%d.mat\n'],detno);
%%
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
