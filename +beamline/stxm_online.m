% Call function without arguments for instructions on how to use it

% Filename: $RCSfile: stxm_online.m,v $
%
% $Revision: 1.16 $  $Date: 2011/04/04 17:03:48 $
% $Author:  $
% $Tag: $
%
% Description:
% plot a STXM scan
%
% Note:
% Call without arguments for a brief help text.
%
% Dependencies: 
% - image_read
%
% history:
%
% April 4th 2011:
% do not normalize the dark field since this is problematic for SAXS with a
% beam stop
%
% September 29th 2010:
% include changes by Martin Dierolf and Joan Vila in the standard version
% of stxm_online
%
% December 10th 2008:
% add bug-fixes and suggestions from Martin Dierolf:
% DirPerLine parameter could not be set via the command line,
% BurstMode flag was always active, is now coupled to dir_per_line, 
% new Parameter ZeroOrderR
%
% September 5th 2008:
% use compile_x12sa_filename, 
% plot as 2x2 sub figures
%
% June 14th 2008: 1st documented version based on work 

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

function [varargout] = stxm_online(first_scan_number, Ny, varargin)
import beamline.pilatus_valid_pixel_roi
import io.image_read
import plotting.image_show
import utils.compile_x12sa_filename
import utils.find_files

% set default values
% Pilatus 2M
detector_number = 1;
% single directory or directory per line format
dir_per_line = 1;
% figure number for display
fig_no = 2;
% number of points along a scan line, 0 for automatic determination from
% the first line
Nx = 0;
% size of the regio of interest
roi_dim = 128;
% automatic determination of the center position 
cen_x = 0;
cen_y = 0;
% dark field integration starting radius
dark_field_r = 20;
% radius of excluded area around center
zero_order_r = 0;
% calculate the first moment rather than a Fourier transform to get the
% differential phase contrast
first_moment = 1;
% use additionally differentiation of the integrated phase
integrated_phase = 1;
% do not update the plot every line to save some time
update_interval = 3;
% save resulting figure
figure_dir = '~/Data10/analysis/online/stxm/figures/';
% save the resulting data
data_dir = '~/Data10/analysis/online/stxm/data/';
% valid pixel mask
filename_valid_mask = '~/Data10/analysis/data/pilatus_valid_mask.mat';

phase = [];
gx = [];
gy = [];

full_screen_position_integrated_phase = [ 5         525        1201         420];
print_a4_position_integrated_phase =  [ 5   525   743   420 ];
full_screen_position_standard = [ 5         109        1201         836];
print_a4_position_standard = [ 5   109   743   836 ];

% check minimum number of input arguments
if (nargin < 2)
    fprintf('Usage:\n')
    fprintf('[trans,dpcx,dpcy,df]=%s(<(first) scan number>, <no. of scan lines> [[,<name>,<value>] ...]);\n',...
            mfilename);
    fprintf('The optional <name>,<value> pairs are:\n');
    fprintf('''DetectorNumber'',<1-Pilatus 2M, 2-Pilatus 300k, 3-Pilatus 100k>\n');
    fprintf('''Nx'',<no. of points per line>         default is %d (0 means automatic determination from first scan line)\n',Nx);
    fprintf('''ROIdim'',<no. of points>              region of interest used for data analysis, default is %d\n',roi_dim);
    fprintf('''CenX'',<point>                        0 means automatic determination, default is %d\n',cen_x);
    fprintf('''CenY'',<point>                        0 means automatic determination, default is %d\n',cen_y);
    fprintf('''DarkFieldR'',<min. radius>            dark field integration starts at this radius, default is %.0f\n',dark_field_r);
    fprintf('''FigNo'',<integer value>               figure number for data display, default is %d\n',fig_no);
    fprintf('''DirPerLine'',<0-no,1-yes>             separate directory for each scan line, default is %d\n',dir_per_line);
    fprintf('''ZeroOrderR'', <min. radius>           pixel values inside this radius are set to zero, default is %d\n', zero_order_r);
    fprintf('''FirstMoment'',<0-no,1-yes>            calculate the first moment rather than a Fourier transform to get the differential phase contrast, default is %d\n',first_moment);
    fprintf('''IntegratedPhase'',<0-no,1-yes>        differentiate additionally the sum signal and re-differentiate it, default is %d\n',integrated_phase);
    fprintf('''UpdateInterval'',<integer N>          update the plot each Nth line, default is %d\n',update_interval);
    fprintf('''FigureDir'',''directory''               save the resulting plot in eps, jpeg and Matlab fig format, '''' for no saving, default is %s\n',figure_dir);
    fprintf('''DataDir'',''directory''                 save the resulting data as Matlab file, '''' for no saving, default is %s\n',data_dir);
    fprintf('''FilenameValidMask'',<path and filename> Matlab file with the valid pixel indices ind_valid, [] for no valid pixel mask,\n');
    fprintf('                                          default is %s\n',filename_valid_mask);
    fprintf('Additional <name>,<value> pairs recognized by compile_x12sa_filename and by image_read can be specified. Please call them for an overview\n');
    fprintf('\n');
    error('At least the (first) scan number and the number of scan lines have to be specified as input parameter.');
end

% accept cell array with name/value pairs as well
no_of_in_arg = nargin;
if (nargin == 3)
    if (isempty(varargin))
        % ignore empty cell array
        no_of_in_arg = no_of_in_arg -1;
    else
        if (iscell(varargin{1}))
            % use a filled one given as first and only variable parameter
            varargin = varargin{1};
            no_of_in_arg = 2 + length(varargin);
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
        case 'DetectorNumber'
            detector_number = value;
        case 'Nx' 
            Nx = value;
        case 'ROIdim' 
            roi_dim = value;
        case 'CenX' 
            cen_x = value;
        case 'CenY' 
            cen_y = value;
        case 'DarkFieldR' 
            dark_field_r = value;
        case 'FigNo' 
            fig_no = value;
        case 'DirPerLine'
            dir_per_line = value;
        case 'ZeroOrderR'
            zero_order_r = value;
        case 'FirstMoment'
            first_moment = value;
        case 'IntegratedPhase'
            integrated_phase = value;
        case 'UpdateInterval'
            update_interval = value;
        case 'FilenameValidMask' 
            filename_valid_mask = value;
        otherwise
            vararg{end+1} = name; %#ok<AGROW>
            vararg{end+1} = value; %#ok<AGROW>
    end
end

% pass some parameters to image_show
vararg(11:(end+10)) = vararg;
vararg{ 1} = 'RetryReadSleep';
vararg{ 2} = 5.0;
vararg{ 3} = 'RetryReadMax';
vararg{ 4} = 5;
vararg{ 5} = 'ErrorIfNotFound';
vararg{ 6} = 0;
% vararg{ 7} = 'BurstMode';
% if (dir_per_line)
%     vararg{ 8} = 1;
% else
%     vararg{ 8} = 0;
% end
vararg{7} = 'UnhandledParError';
vararg{8} = 0;
vararg{9} = 'DetectorNumber';
vararg{10} = detector_number;

% region of interest index in each dimension
roi_rel_ind = -round(0.5*roi_dim):(round(0.5*roi_dim)-1);

% load the indices of valid pixels
if ((~isempty(filename_valid_mask)) && (exist(filename_valid_mask,'file')))
    fprintf('loading the valid pixel mask %s\n',filename_valid_mask);
    load(filename_valid_mask);
end

% wait for the data to be available
scan_no_check = first_scan_number;
if ((dir_per_line) && (Ny > 1))
    scan_no_check = scan_no_check +1;
end
filename_mask = compile_x12sa_filename(scan_no_check,0,'DetectorNumber',detector_number);
[~, fnames] = find_files(filename_mask);
data_available = (~isempty(fnames));
if (~data_available)
    fprintf('Waiting for %s to become available.\n',filename_mask);
    while (~data_available);
        pause(1);
        [~, fnames] = find_files(filename_mask);
        data_available = (~isempty(fnames));
    end
end

% check that number of points per line determination will be possible

% determine number of points per line
if (Nx <= 0)
    if (~dir_per_line)
        error('The number of points per line can only automatically be determined if separate scan directories are used for each line.');
    end
    vararg_remain = vararg;
    vararg_remain(3:(end+2)) = vararg_remain;
    vararg_remain{1} = 'SubExpWildcard';
    vararg_remain{2} = 1;
    [fmask,vararg_remain] = ...
        compile_x12sa_filename(first_scan_number,0,vararg_remain); %#ok<NASGU>
    Nx = length(dir(fmask));
    if (Nx < 1)
        error('No matching files found for %s',fmask);
    end
end
fprintf('%d lines with %d points per line in\n',Ny,Nx);


if (integrated_phase)
    figure(fig_no +1);
    hold off;
    clf;
    % print as layed out on the screen, i.e., preserve aspect ratio
    set(gcf,'PaperPositionMode','auto');
    % paper size
    set(gcf,'PaperType','A4');
    % background color
    set(gcf,'Color','white');
    % resize and position
    set(gcf,'Position',full_screen_position_integrated_phase);

    colormap(bone(256));
end


figure(fig_no);
hold off;
clf;
% print as layed out on the screen, i.e., preserve aspect ratio
set(gcf,'PaperPositionMode','auto');
% paper size
set(gcf,'PaperType','A4');
% background color
set(gcf,'Color','white');
% resize and position
set(gcf,'Position',full_screen_position_standard);

colormap(bone(256));


% STXM display loop
point_no = 0;
scan_number = first_scan_number;

frame = [];
for ii=Ny:-1:1
    sub_exp_no = 0;
    for jj=Nx:-1:1
        if (dir_per_line)
            vararg_remain = vararg;
            vararg_remain(3:(end+2)) = vararg_remain;
            vararg_remain{1} = 'SubExpNo';
            vararg_remain{2} = sub_exp_no;
            [filename,vararg_remain] = ...
                compile_x12sa_filename(scan_number,0,vararg_remain); 
        else
            [filename,vararg_remain] = ...
                compile_x12sa_filename(scan_number,point_no,vararg); 
        end
        last_frame = frame;
        [frame,vararg_remain] = image_read(filename,vararg_remain);
        if (isempty(frame.data))
            fprintf('%s not found, repeating the previous frame\n',filename);
            frame = last_frame;
        end
        if (~isempty(vararg_remain))
            vararg_remain
            error('There are unhandled parameters.');
        end
        
        if (point_no == 0)            
            trans = zeros(Ny,Nx);
            dpcx = trans;
            dpcy = trans;
            df = trans;
      
            if ((cen_x <= 0) || (cen_y <= 0))
                [cx, cy] = find_center(frame.data);
                fprintf('Beam cemter guess (x,y) = (%d,%d)\n',cx,cy);
                if (cen_x <= 0)
                    cen_x = cx;
                end
                if (cen_y <= 0)
                    cen_y = cy;
                end
            end
            
            roi_x_ind = cen_x + roi_rel_ind;
            if ((roi_x_ind(1) < 1) || (roi_x_ind(end) > size(frame.data,2)))
                error('Region of interest out of range in x\n');
            end
            roi_y_ind = cen_y + roi_rel_ind;
            if ((roi_y_ind(1) < 1) || (roi_y_ind(end) > size(frame.data,1)))
                error('Region of interest out of range in y\n');
            end

            [yy,xx] = meshgrid(roi_rel_ind,roi_rel_ind);
            [~, rho] = cart2pol(xx,yy);

            ind_df = find((rho > dark_field_r) & (rho < roi_rel_ind(end)));
            
            if (~isempty(filename_valid_mask))
                % in case of less than full detector readout cut out the right part of 
                % the valid pixel mask
                valid_mask = pilatus_valid_pixel_roi(valid_mask,'RoiSize',size(frame.data));
            else
                % if the valid pixel mask is not used specify all pixels to
                % be valid
                valid_mask.indices = 1:(size(frame.data,1)*size(frame.data,2));
            end
            
            % calculate the indices of the valid and invalid pixels within
            % the region of interest
            frame_valid = zeros(size(frame.data));
            frame_valid(valid_mask.indices) = 1;
            frame_valid = frame_valid(roi_y_ind,roi_x_ind);
            ind_invalid = find(frame_valid == 0);
%             ind_valid = find(frame_valid ~= 0);
            ind_df = setdiff(ind_df,ind_invalid);           
        end

        % cut out the region of interest
        frame_roi = frame.data(roi_y_ind,roi_x_ind);
        frame_roi(ind_invalid) = 0;
        
        % set central part of detector frame to zero, if specified
        if (zero_order_r> 0)
            frame_roi(rho<zero_order_r) = 0; %min(frame_roi(:));
        end

        % data analysis for the current point
        if (first_moment)
            [tr,px,py] = stxm_pt2(frame_roi);
        else
            [tr,px,py] = stxm_pt(frame_roi);
        end
        trans(ii,jj) = tr;
        dpcx(ii,jj) = px;
        dpcy(ii,jj) = py;
      
%         df(ii,jj) = sum(frame_roi(ind_df)) / sum(sum(frame_roi(ind_valid)));
        df(ii,jj) = sum(frame_roi(ind_df));
      
        point_no = point_no +1;
        sub_exp_no = sub_exp_no +1;
    end
   
    % plot linewise each update_interval-th line
    if (Nx > 1) && (Ny > 1)
        if ((ii == Ny) || (rem(ii,update_interval) == 1) || (ii == 1))
            if(gcf ~= fig_no)
                figure(fig_no);
            end
            iv = 2;
            ih = 2;
            colormap(bone(256));

            subplot(iv,ih,1); 
            imagesc(trans); 
            axis xy; axis equal; axis tight;
            colorbar; 
            axis_min = min(min(trans(trans ~= 0)));
            if (isnan(axis_min))
                axis_min = 0;
            end
            axis_max = max(max(trans(trans ~= 0)));
            if (isnan(axis_max))
                axis_max = 0;
            end            
            caxis([(axis_min-.0001) (axis_max+.0001)]);
            title_str = [ 'transmission #' num2str(first_scan_number,'%d') ];
            if (dir_per_line)
               title_str = [ title_str '-' num2str(first_scan_number+Ny-1,'%d') ]; %#ok<AGROW>
            end
            title_str = sprintf('%s (detector %d)',title_str,detector_number);
            title(title_str);

            subplot(iv,ih,2); 
            imagesc(df); 
            axis xy; axis equal; axis tight;
            colorbar; 
            axis_min = min(min(df(df ~= 0)));
            if (isnan(axis_min))
                axis_min = 0;
            end
            axis_max = max(max(df(df~=0)));
            if (isnan(axis_max))
                axis_max = 0;
            end            
            caxis([(axis_min-.0001) (axis_max+.0001)]);
            title('dark field');

            subplot(iv,ih,3); 
            imagesc(dpcx); 
            axis xy; axis equal; axis tight;
            colorbar; 
            axis_min = min(min(dpcx(dpcx ~= 0)));
            if (isnan(axis_min))
                axis_min = 0;
            end
            axis_max = max(max(dpcx(dpcx~=0)));
            if (isnan(axis_max))
                axis_max = 0;
            end            
            caxis([(axis_min-.0001) (axis_max+.0001)]);
            title('DPC x');

            subplot(iv,ih,4); 
            imagesc(dpcy); 
            axis xy; axis equal; axis tight;
            colorbar; 
            axis_min = min(min(dpcy(dpcy ~= 0)));
            if (isnan(axis_min))
                axis_min = 0;
            end
            axis_max = max(max(dpcy(dpcy~=0)));
            if (isnan(axis_max))
                axis_max = 0;
            end            
            caxis([(axis_min-.0001) (axis_max+.0001)]);
            title('DPC y');

            drawnow;
        end
    end
   
    if (dir_per_line)
       scan_number = scan_number +1;
    end
end

% store return arguments
if nargout > 0
  varargout{1} = trans;
end
if nargout > 1
  varargout{2} = dpcx;
end
if nargout > 2
  varargout{3} = dpcy;
end
if nargout > 3
    varargout{4} = df;
end

if nargout > 4
    varargout{5} = phase;
end

if nargout > 5
    varargout{6} = gx;
end

if nargout > 6
    varargout{7} = gy;
end

if (integrated_phase)
    % calculate the integrated phase from the differential phase contrast
    % in horizontal and vertical direction
    phase = phase_from_dpc(dpcx,dpcy, 'fourier');
    
    % calculate the 1D differential phase contrast from the integrated
    % phase
    [gx, gy] = gradient(phase);
    
    figure(fig_no +1);
    iv = 1;
    ih = 3;
    colormap(bone(256));

    subplot(iv,ih,1); 
    imagesc(phase); 
    axis xy; axis equal; axis tight;
    colorbar; 
    axis_min = min(phase(phase ~= 0));
    if (isnan(axis_min))
        axis_min = 0;
    end
    axis_max = max(phase(phase~=0));
    if (isnan(axis_max))
        axis_max = 0;
    end    
    caxis([(axis_min-.0001) (axis_max+.0001)]);
    title_str = [ 'integrated phase #' num2str(first_scan_number,'%d') ];
    if (dir_per_line)
       title_str = [ title_str '-' num2str(first_scan_number+Ny-1,'%d') ]; 
    end
    title_str = sprintf('%s (detector %d)',title_str,detector_number);
    title(title_str);

    subplot(iv,ih,2); 
    imagesc(gx); 
    axis xy; axis equal; axis tight;
    colorbar; 
    axis_min = min(gx(gx ~= 0));
    if (isnan(axis_min))
        axis_min = 0;
    end
    axis_max = max(gx(gx ~= 0));
    if (isnan(axis_max))
        axis_max = 0;
    end    
    caxis([(axis_min-.0001) (axis_max+.0001)]);
    title('DPC x from integrated phase');

    subplot(iv,ih,3); 
    imagesc(gy);
    axis xy; axis equal; axis tight;
    colorbar; 
    axis_min = min(gy(gy ~= 0));
    if (isnan(axis_min))
        axis_min = 0;
    end
    axis_max = max(gy(gy ~= 0));
    if (isnan(axis_max))
        axis_max = 0;
    end    
    caxis([(axis_min-.0001) (axis_max+.0001)]);
    title('DPC y from integrated phase');

    drawnow;

end


% file name for saving
filename = sprintf('stxm_scans_%d_%05d-%05d',detector_number,...
    first_scan_number,first_scan_number+Ny-1);


% save figures
if (~isempty(figure_dir))
    figure(fig_no);

    % create output directories and write the plot in different formats
    if (~exist(figure_dir,'dir'))
        mkdir(figure_dir)
    end
    if ((figure_dir(end) ~= '/') && (figure_dir(end) ~= '\')) 
        figure_dir = [ figure_dir '/' ];
    end
    fprintf('output directory for figures is %s\n',figure_dir);
    
    % resize to a smaller width as print layout
    set(gcf,'Position',print_a4_position_standard);
    
    subdir = [ figure_dir 'jpg/' ];
    if (~exist(subdir,'dir'))
        mkdir(subdir);
    end
    fprintf('saving %s.jpg\n',filename);
    print('-djpeg','-r300',[subdir filename '.jpg'] );
    
    subdir = [ figure_dir 'eps/' ];
    if (~exist(subdir,'dir'))
        mkdir(subdir);
    end
    fprintf('saving %s.eps\n',filename);
    print('-depsc','-r1200',[subdir filename '.eps'] );

    % resize to full screen
    set(gcf,'Position',full_screen_position_standard);

    subdir = [ figure_dir 'fig/' ];
    if (~exist(subdir,'dir'))
        mkdir(subdir);
    end
    fprintf('saving %s.fig\n',filename);
    hgsave([subdir filename '.fig']);    

    if (integrated_phase)
        figure(fig_no +1);

        subdir = [ figure_dir 'jpg/' ];
        if (~exist(subdir,'dir'))
            mkdir(subdir);
        end
        
        % resize to a smaller width as print layout
        set(gcf,'Position',print_a4_position_integrated_phase);
        
        fprintf('saving %s_integrated_phase.jpg\n',filename);
        print('-djpeg','-r300',[subdir filename '_integrated_phase.jpg'] );

        subdir = [ figure_dir 'eps/' ];
        if (~exist(subdir,'dir'))
            mkdir(subdir);
        end
        fprintf('saving %s_integrated_phase.eps\n',filename);
        print('-depsc','-r1200',[subdir filename '_integrated_phase.eps'] );

        % resize to a smaller width as print layout
        set(gcf,'Position',full_screen_position_integrated_phase);
                
        subdir = [ figure_dir 'fig/' ];
        if (~exist(subdir,'dir'))
            mkdir(subdir);
        end
        fprintf('saving %s_integrated_phase.fig\n',filename);
        hgsave([subdir filename '_integrated_phase.fig']);    
    end
end


% save resulting data
if (~isempty(data_dir))
    if ((data_dir(end) ~= '/') && (data_dir(end) ~= '\')) 
        data_dir = [ data_dir '/' ];
    end
    
    % create output directory
    if (~exist(data_dir,'dir'))
        mkdir(data_dir)
    end
    
    % save data
    fprintf('saving %s.mat\n',[data_dir filename]);
    if (integrated_phase)
        save([data_dir filename],'trans','dpcx','dpcy','df', 'phase', 'gx', 'gy');
    else
        save([data_dir filename],'trans','dpcx','dpcy','df');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cx, cy] = find_center(f)

f = medfilt2(f,[5 5]);

[~, cx] = max(sum(f,1));
[~, cy] = max(sum(f,2));

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tr, px, py] = stxm_pt(a)

    persistent c1 c2 s1 s2 sz

    if (isempty(sz)) || (any(sz ~= size(a)))
        sz = size(a);
        c1 = -cos(2*pi*(0:sz(1)-1)/sz(1));
        s1 =  sin(2*pi*(0:sz(1)-1)/sz(1));
        c2 = -cos(2*pi*(0:sz(2)-1)/sz(2));
        s2 =  sin(2*pi*(0:sz(2)-1)/sz(2));
    end

    a1 = sum(a,1);
    a2 = sum(a,2)';

    tr = sum(a1);
    px = atan2(sum(a1.*c1), sum(a1.*s1));
    py = atan2(sum(a2.*c2), sum(a2.*s2)); 

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tr, px, py ] = stxm_pt2(a)

    persistent x y sz

    if (isempty(sz)) || (any(sz ~= size(a)))
        sz = size(a);
        % masking out the invalid pixels is done by setting the
        % corresponding intensities to zero before calling this function
        [y,x] = ndgrid((0:sz(1)-1)-sz(1)/2, (0:sz(1)-1)-sz(1)/2);
%         x2 = x.^2;
    end

    tr = sum(sum(a));
    px = sum(sum(a.*x))/tr;
    py = sum(sum(a.*y))/tr;
%     p2 = (sum(a1.*x2)/tr + sum(a2.*x2)/tr - px^2 - py^2);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = phase_from_dpc(dpcx,dpcy,varargin)
%
% Integrates the phase from a combination of x and y gradients.
% phase_from_dpc(dpcx,dpcy,'fourier') uses the Fourier method (default),
% phase_from_dpc(dpcx,dpcy,'finitdiff') uses a finite difference method.

if nargin > 2
    method = varargin{1};
else
    %method = 'fourier';
    method = 'finitediff';
end

px = -dpcy;
py = -dpcx;

sz = size(px);

switch lower(method)
    case 'fourier'
        f = zeros(2*sz);
        f(1:sz(1),1:sz(2)) = px + 1i*py;
        f(1:sz(1),sz(2)+1:end) = fliplr(px + 1i*py);
        f(sz(1)+1:end,1:sz(2)) = flipud(px + 1i*py);
        f(sz(1)+1:end,sz(2)+1:end) = rot90(px + 1i*py,2);
        [x1,x2] = ndgrid(-sz(1):(sz(1)-1),-sz(2):(sz(2)-1));
        q1 = pi*fftshift(x1)/sz(1);
        q2 = pi*fftshift(x2)/sz(2);
        qc = q2 - 1i*q1;
        inv_qc = 1./qc;
        inv_qc(1,1) = 0;
        nf = ifftn(fftn(f).*inv_qc);
        p = real(nf(1:sz(1),1:sz(2)));
    case 'finitediff'
        ggx = pgradient(dpcx);
        [~, ggy] = pgradient(dpcy);
        f = .25*(ggx + ggy);
        ta = zeros(sz);
        for i = 1:10000
             ta = ta + (pdel2(ta) - f);

             % Zero boundary conditions
             %ta(1,:) = 0;
             %ta(:,1) = 0;
             %ta(end,:) = 0;
             %ta(:,end) = 0;

             % Zero normal gradient boundary condition
             ta(1,:) = ta(2,:);
             ta(:,1) = ta(:,2);
             ta(end,:) = ta(end-1,:);
             ta(:,end) = ta(:,end-1);
    
             if mod(i,1000)==0
                 figure(1); imagesc(real(ta)); colormap(bone(256)); colorbar; drawnow;
             end

             p = ta;
        end
end

