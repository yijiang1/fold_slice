% Call function without arguments for instructions on how to use it

% Filename: $RCSfile: prep_integ_masks.m,v $
%
% $Revision: 1.9 $  $Date: 2016/01/21 14:51:57 $
% $Author:  $
% $Tag: $
%
% Description:
% prepare masks for the radial integration of SAXS patterns
%
% Note:
% Call without arguments for a brief help text.
%
% Dependencies:
% - image_read
% - pilatus_valid_pixel_roi
%
% history:
%
% September 4th 2009:
% correct in help text one of the RadiusFrom to RadiusTo
%
% May 9th 2008: 1st documented version
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

function [ integ_masks ] = prep_integ_masks(filename, center_xy, varargin)
import beamline.pilatus_valid_pixel_roi
import io.image_read
import plotting.display_valid_mask
import utils.pixel_to_q

% set number of radii
no_of_radii = 0;
% number of angular segments per radius
no_of_segments = 1;
% pixel size in mm
pixel_size_mm = [];%.172;
% detector distance in mm
det_dist_mm = [];%2000;
calculate_q=0; %only calculate q if exact detector distance is given
% wavelength (unit inconsequential, will be reflected in q)
lambda = [];%1;
% output directory for the masks
out_dir = '~/Data10/analysis/data/';
filename_valid_mask = [ out_dir 'pilatus_valid_mask.mat' ];
filename_integ_masks = [ out_dir 'pilatus_integration_masks.mat' ];
% output figure number
fig_no = 240;
% save integration masks
save_data = 1;
% display valid pixel mask
display_valid_mask_flag = 1;
% detector number
det_no= 1;
% angular range to be excluded to cut out the beam stop
bs_angle_from = 0;
bs_angle_to = 0;

% check minimum number of input arguments
if (nargin < 2)
    fprintf('\nUsage:\n');
    fprintf('[integ_masks]=%s( filename, center_xy [[,<name>,<value>]...]);\n',mfilename);
    fprintf('Prepare the masks for an efficient radial integration.\n');
    fprintf('The optional angular range in degree can be used to cut out a beam stop.\n');
    fprintf('Angle 0 is horizontally to the left, positive in counterclockwise direction.\n');
    fprintf('The specified data file is loaded and some of the integration masks are plotted into that frame.\n');
    fprintf('\n');
    fprintf('The optional <name>,<value> pairs are:\n');
    
    fprintf('''NormalXY'',[x y]                    pixel coordinates, from where the detector normal points\n');
    fprintf('                                    toward the sample.  Default is equal to center_xy\n');
    fprintf('''PixelSize_mm'',<double>             pixel size in mm, default is %.3f\n',pixel_size_mm);
    fprintf('''DetDist_mm'',<double>               detector distance in mm, default is %.1f\n',det_dist_mm);
    fprintf('''Wavelength'',<double>               wavelength. The units chosen here will determine the units of q\n');
    fprintf('                                    The defaults is %.1f\n',lambda);
    fprintf('''NoOfRadii'',<integer>               radial integration start radius, default is %d\n',no_of_radii);
    fprintf('        or ,<vector>,              defining the limits of radius bins\n');
    fprintf('''NoOfSegments'',<integer>            Number of angular segments. If an integer, this number of equally wide azimuthal\n')
    fprintf('                                    bins over 360 degrees are created. default is %d\n',no_of_segments);
    fprintf('        or ,<vector>,              defining the limits of angular bins\n');
    fprintf('''SaveData'',<0-no,1-yes>             save the integration masks, default is %d\n',save_data);
    fprintf('''FilenameValidMask'',<path and filename>  Matlab file with the valid pixel indices ind_valid,\n');
    fprintf('                                    default is %s\n',filename_valid_mask);
    fprintf('''FilenameIntegMasks'',<path and filename> output file name for the structure integ_masks,\n');
    fprintf('                                    default is %s\n',filename_integ_masks);
    fprintf('''FigNo'',<integer>                   number of the figure in which the result is displayed\n');
    fprintf('''DetNo'',<integer>                   number of detector 1 for SAXS and 2 for WAXS\n');
    fprintf('                                     Default is 1 (SAXS)\n');
    fprintf('''BeamstopAngleFrom'',<float>         exclude an angular region from the integration, default for the start value is %d\n',...
        bs_angle_from);
    fprintf('''BeamstopAngleTo'',<float>           exclude an angular region from the integration, default for the end value is %d\n',...
        bs_angle_to);
    fprintf('\n');
    fprintf('\n');
    fprintf('The file name should be the name of a single file without wildcards\n');
    fprintf('that is displayed as an example.\n');
    fprintf('The image file has no other function beyond being displayed as example.\n');
    fprintf('Example:\n');
    fprintf('[integ_masks]=%s(''~/Data10/pilatus/image.cbf'',[512 512]);\n',...
        mfilename);
    
    error('At least the filename and the beam center have to be specified as input parameter.');
end

% check number of center coordinates
if (length(center_xy) ~= 2)
    error('The beam center needs to be specified as a two component vector [cen_x cen_y].\n');
end
center_x = center_xy(1);
center_y = center_xy(2);
norm_x = center_x;
norm_y = center_y;

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
        case 'NormalXY'
            if (numel(value)==2)
                norm_x = value(1);
                norm_y = value(2);
            end
        case 'PixelSize_mm'
            pixel_size_mm = value;
        case 'DetDist_mm'
            det_dist_mm = value;
            calculate_q=1;
        case 'Wavelength_nm'
            lambda = value;
        case 'Wavelength'
            lambda = value/10;
        case 'NoOfRadii'
            no_of_radii = value;
        case 'NoOfSegments'
            no_of_segments = value;
        case 'FilenameValidMask'
            filename_valid_mask = value;
        case 'FilenameIntegMasks'
            filename_integ_masks = value;
        case 'SaveData'
            save_data = value;
        case 'DisplayValidMask'
            display_valid_mask_flag = value;
        case 'FigNo'
            fig_no = value;
        case 'DetNo'
            det_no = value;
        case 'BeamstopAngleFrom' 
            bs_angle_from = value;
        case 'BeamstopAngleTo' 
            bs_angle_to = value; 
        otherwise
            vararg{end+1} = name;
            vararg{end+1} = value;
    end
end



% MGS - The fact that the detector number is dictating whether to have a
% radius variable or q is not ideal. Additional arguments should be given
% for this such as calculate_q or save_radius_var.
% initialize return arguments
if (det_no == 1)||(det_no == 3)
    calculate_radius = true;
elseif (det_no == 2)
    calculate_radius = false;
else
    error('Only det_no 1, 2, and 3 are recognized')
end
if ( calculate_radius ) && (calculate_q)
integ_masks = struct('radius',[], 'indices',[], 'norm_sum', [],'q',[]);
elseif (calculate_radius) && (calculate_q == 0)
  integ_masks = struct('radius',[], 'indices',[], 'norm_sum', []); 
elseif (~calculate_radius)
  integ_masks = struct('indices',[], 'norm_sum', [],'q',[]); 
end


% check radius
%if (exist('r_from','var'))
if (size(no_of_radii)>1)
    if (no_of_radii(1) < 1)
        %if (r_from < 1)
        error('The minimum radius is 1, %d is invalid.',r_from);
    end
    %end
    %if (exist('r_from','var') && exist('r_to','var'))
    if  (no_of_radii(end) < no_of_radii(1))
        %if ((r_to ~= 0) && (r_to < r_from))
        error('The maximum radius must be greater than the minimum one, %d is invalid.\n',no_of_radii(end));
    end
end

% check number of angular segments
if (no_of_segments < 1)
    error('The number of angular segments must be at least 1');
end
if (numel(no_of_segments)>1)
    angular_segments = no_of_segments;
    no_of_segments = numel(no_of_segments)-1;
else
    angular_segments = 360/no_of_segments * (0:no_of_segments);
end
angular_segments = mod(angular_segments, 360);

% check beamstop region
if ((bs_angle_from < 0.0) || (bs_angle_to > 360.0))
    error('The angular range for the beam stop region is 0 to 360 degree.');
end
if (bs_angle_to < bs_angle_from)
    error('The maximum beam stop angle must be less than or equal to the minimum one.\n');
end

% load the indices of valid pixels
fprintf('loading the valid pixel mask %s\n',filename_valid_mask);
load(filename_valid_mask);
dim_x = valid_mask.framesize(2);
dim_y = valid_mask.framesize(1);


if (~isempty(filename))
    % load test frame
    frame = image_read(filename,vararg);
    % select the first frame for display
    frame.data = frame.data(:,:,1);
    % in case of less than full detector readout cut out the right part of
    % the valid pixel mask
    valid_mask = pilatus_valid_pixel_roi(valid_mask,'RoiSize',size(frame.data));
    dim_x = size(frame.data,2);
    dim_y = size(frame.data,1);
end

% plot valid pixel mask
if (display_valid_mask_flag)
    figure(fig_no);
    vpm = zeros(dim_y,dim_x);
    vpm(valid_mask.indices) = 1;
    imagesc(vpm);
    axis xy;
    axis equal;
    axis tight;
    title('valid pixels');
    set(gcf,'Name','valid pixels');
    drawnow;
end

if calculate_radius
    %if (exist('r_to','var'))
    % choose maximum radius, if specified via r_to=0
    if (no_of_radii < 1)
        no_of_radii = max( [ sqrt(center_x^2+center_y^2) ...
            sqrt((dim_x-center_x)^2+center_y^2) ...
            sqrt(center_x^2+(dim_y-center_y)^2) ...
            sqrt((dim_x-center_x)^2+(dim_y-center_y)^2) ] );
    end
    if size(no_of_radii)==1
        no_of_radii=1:1:no_of_radii;
    end
end

fprintf('preparing the integration masks ...\n');

% create an array of the (x,y) coordinates relative to the beam center and
% convert it to polar coordinates
%
if calculate_radius % For SAXS detector - MGS, should be fixed, why is it neded different calculation for different detectors?
    % angular range to be excluded to cut out the beam stop

    [ x, y ] = meshgrid( (1:dim_x)-center_x, (1:dim_y)-center_y );
    [ theta, rho ] = cart2pol( x, y );
    % convert angular range from -pi/pi to 0/360
    theta = (theta/pi +1) * 180.0;
    
    % prepare circular masks of the integer width r_step (in pixel)
    integ_masks.radius = no_of_radii;
    r_step=no_of_radii(2)-no_of_radii(1);
    if calculate_q
        integ_masks.q = pixel_to_q(no_of_radii,pixel_size_mm,det_dist_mm, 12.39852/lambda);
    end
    no_of_radii = length(integ_masks.radius);
    integ_masks.indices = cell( no_of_radii, no_of_segments );
    integ_masks.norm_sum = zeros( no_of_radii, no_of_segments );
    seg_inds = cell(no_of_segments,1);
    for ind_seg = 1:no_of_segments
        seg_from = angular_segments(ind_seg);
        seg_to   = angular_segments(ind_seg+1);
        if (seg_from >= seg_to)
            ind_curr = find( ((theta > seg_from) | (theta <= seg_to) ) & ...
                ((theta <= bs_angle_from) | (theta >= bs_angle_to)) );
        else
            ind_curr = find( ((theta > seg_from) & (theta <= seg_to) ) & ...
                ((theta <= bs_angle_from) | (theta >= bs_angle_to)) );
        end
        % only take valid pixels into account
        ind_curr = intersect(ind_curr, valid_mask.indices);
        seg_inds{ind_seg} = ind_curr;
    end
    
    for ind_r=1:no_of_radii
        if (rem(ind_r,100) == 0)
            fprintf('%4d / %d',ind_r,no_of_radii);
            if (ind_r <= no_of_radii-100)
                fprintf(', ');
            end
        end
        r_inds = find( (rho >= integ_masks.radius(ind_r)) & ...
            (rho < integ_masks.radius(ind_r)+r_step) );
        for ind_seg = 1:no_of_segments
            integ_masks.indices{ind_r, ind_seg} = intersect( r_inds, seg_inds{ind_seg} );
            % calculate the normalization value (sum of the pixels within the mask)
            integ_masks.norm_sum(ind_r, ind_seg) = ...
                length( integ_masks.indices{ind_r, ind_seg} );
        end
    end
    fprintf('\n');
    
else 
    [ x, y ] = meshgrid( (1:dim_x)-norm_x, (1:dim_y)-norm_y );
    if (norm_x == center_x && norm_y == center_y)
        [ theta, rho ] = cart2pol( x, y );
        q = 4*pi/lambda*sin(atan2(rho,det_dist_mm/pixel_size_mm)/2);
        % convert angular range from -pi/pi to 0/360
        theta = theta/pi*180.0;
    else
        if (norm_x ~= center_x)
            angle = atan((norm_x - center_x) / (det_dist_mm/pixel_size_mm));
            z = -x*sin(angle) + det_dist_mm/pixel_size_mm*cos(angle);
            x =  x*cos(angle) + det_dist_mm/pixel_size_mm*sin(angle);
        else
            fprintf('not implemented yet!!!\n');
            exit
        end
        q = 4*pi/lambda*sin(atan2(sqrt(x.^2 + y.^2),z)/2);
        theta = atan2(y,x)/pi*180;
    end
    t_1d = reshape(theta(valid_mask.indices),1,[]);
    q_1d = reshape(q(valid_mask.indices),1,[]);
    
    t_ed = linspace(     -180,      180,1e0+1);
    integ_masks.theta = t_ed(1:end-1);
    integ_masks.theta_end = t_ed(end);
    q_ed = linspace(min(q_1d),max(q_1d),1e3+1);
    integ_masks.q = q_ed(1:end-1);
    integ_masks.q_end = q_ed(end);
    
    [~,t_bin] = histc(t_1d,t_ed);
    [~,q_bin] = histc(q_1d,q_ed);
    
    integ_masks.indices = cell(numel(q_ed)-1,numel(t_ed)-1);
    integ_masks.norm_sum = zeros(size(integ_masks.indices));
    for q_i=1:numel(q_ed)-1
        for t_i=1:numel(t_ed)-1
            integ_masks.indices{q_i,t_i} = ...
                valid_mask.indices(and(q_bin==q_i,t_bin==t_i));
            integ_masks.norm_sum(q_i,t_i) = numel(integ_masks.indices{q_i,t_i});
        end
    end
end


% save integration masks
if (save_data)
    fprintf('Saving center_xy, no_of_segments, integ_masks to %s\n',...
        filename_integ_masks);
    if angular_segments(end) == 0
        angular_segments(end) = 360;
    end
    phi_det = (angular_segments(2:end) + angular_segments(1:end-1))/2; %% Center of the angular sector in degrees
    save(filename_integ_masks,'center_xy','no_of_segments','integ_masks','angular_segments','phi_det');
end

% display some integration circles
if (~isempty(filename))
    figure(fig_no+1);
    hold off;
    clf;
    frame_plot = double(frame.data);
    frame_plot( frame_plot < 1 ) = 1;
    
    plot_step = round(length(integ_masks.indices)/50);
    if (plot_step < 2)
        plot_step = 2;
    end
    for (ind_r = 1:plot_step:size(integ_masks.indices,1))
        for (ind_seg = 1:2:no_of_segments)
            frame_plot(integ_masks.indices{ind_r,ind_seg}) = 10^(6*ind_seg/no_of_segments);
        end
    end
    
    imagesc(log10(frame_plot));
    axis xy;
    axis equal;
    axis tight;
    colorbar;
    title([ 'integration segment test plot for ' strrep(filename,'_','\_') ]);
    set(gcf,'Name','integration masks');
end
