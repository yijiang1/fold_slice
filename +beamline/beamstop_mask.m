% Call function without arguments for a detailed explanation of its use

% Filename: $RCSfile: beamstop_mask.m,v $
%
% $Revision: 1.8 $  $Date: 2011/08/23 17:17:53 $
% $Author: $
% $Tag: $
%
% Description:
% remove a polygonic region from the valid pixel mask
%
% Note:
% This is a template. The coordinates of the polygon have to be manually
% edited. 
% Call without arguments for a brief help text.
%
% Dependencies: 
% - image_read
%
% history:
%
% May 19th 2010:
% add XyCoord and xCoord, yCoord command line parameters
%
% May 9th 2008: 1st documented version

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

function [ bmask_ind ] = beamstop_mask(filename,varargin)
import beamline.pilatus_valid_pixel_roi
import beamline.prep_valid_mask
import io.image_read
import plotting.display_valid_mask

% set default values for the variable input arguments:
% valid pixel mask
filename_valid_mask = '~/Data10/analysis/data/pilatus_valid_mask.mat';
% do not update the valid pixel mask
save_data = 0;
% figure number for display
fig_no = 220;
% mask corners
xy_coord = [];  %#ok<NASGU>
x_coord = [];
y_coord = [];

% check minimum number of input arguments
if (nargin < 1)
    display_help(filename_valid_mask,save_data,fig_no);
    error('At least the filename has to be specified as input parameter.');
end

% accept cell array with name/value pairs as well
no_of_in_arg = nargin;
if (nargin == 2)
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
if (rem(no_of_in_arg,2) ~= 1)
    display_help(filename_valid_mask,save_data,fig_no);
    error('The optional parameters have to be specified as ''name'',''value'' pairs');
end


% parse the variable input arguments
vararg_remain = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'SaveData' 
            save_data = value;
        case 'FilenameValidMask' 
            filename_valid_mask = value;
        case 'xyCoord'
            xy_coord = value;
            x_coord = xy_coord(:,1);
            y_coord = xy_coord(:,2);
        case 'xCoord'
            x_coord = value;
        case 'yCoord'
            y_coord = value;
        otherwise
            vararg_remain{end+1} = name; %#ok<AGROW>
            vararg_remain{end+1} = value; %#ok<AGROW>
    end
end

% read file for test display
frame = image_read(filename,vararg_remain);
frame.data = double(frame.data);
dimensions = size(frame.data);
if (numel(dimensions) > 2)
    frame.data = mean(frame.data,3);
    dimensions = size(frame.data);
end

% get indices to pixels within beam stop
if ((isempty(x_coord)) || (isempty(y_coord)))
    bmask_ind = 1:(dimensions(1)*dimensions(2));
else
    [bmask] = uint8(1 - roipoly( dimensions(1), dimensions(2), x_coord, y_coord ));
    bmask_ind = find(bmask == 0);
end

% plot the result
figure(5);
frame_plot = frame.data;
frame_plot(frame_plot < 1) = 1;
% plot the masked region with lower intensity
frame_plot(bmask_ind) = 0.1 * frame_plot(bmask_ind);
imagesc(log10(frame_plot));
axis xy;
axis equal;
axis tight;
colorbar;
title('beamstop mask shape');

% show the current valid pixel mask
display_valid_mask('FilenameValidMask',filename_valid_mask,'FigNo',fig_no+1,...
    'NoHelp',1);
title('current valid pixel mask');

% load ind_valid, the indices of the valid pixels
fprintf('loading %s\n',filename_valid_mask);
load(filename_valid_mask);
% cut out the current region of interest
valid_mask = pilatus_valid_pixel_roi(valid_mask,'RoiSize',size(frame.data));

% remove beam-stop pixels from it
valid_mask.indices = setdiff(valid_mask.indices,bmask_ind); %#ok<NODEF>

if (save_data)
    % create a backup of the mask
    if (exist(filename_valid_mask,'file'))
        filename_valid_mask_backup = [ filename_valid_mask '.bak' ];
        fprintf('Copying the current mask %s to %s\n',filename_valid_mask,...
            filename_valid_mask_backup);
        copyfile(filename_valid_mask,filename_valid_mask_backup);
    end

    % save the updated mask
    fprintf('saving updated mask %s\n',filename_valid_mask);
    save(filename_valid_mask,'valid_mask');

    % display the new mask
    display_valid_mask('FilenameValidMask',filename_valid_mask,'FigNo',fig_no+2,...
    'NoHelp',1);
else
    % mark the valid pixels as 1, leave the invalid at 0
    pframe = zeros(valid_mask.framesize);
    pframe(valid_mask.indices) = 1;

    % plot the result
    figure(fig_no+2);
    imagesc(pframe);
    axis xy;
    axis equal;
    axis tight;
    colorbar;
    title('valid pixels');
    title('updated valid pixel mask (not saved!)');
    set(gcf,'Name','valid pixels');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = display_help(filename_valid_mask,save_data,fig_no)

fprintf('Usage:\n');
fprintf('%s(filename_for_display, [[<name>,<value>],...]);\n',mfilename)
fprintf('The specified file is used to display the beamstop mask with reduced intensity.\n');
fprintf('The optional <name>,<value> pairs are:\n');
fprintf('''xyCoord'',[ x1 y1; x2 y2; ...]      coordinates of the beamstop mask\n');
fprintf('''xCoord'',[ x1 x2 ...]               x-coordinates of the beamstop mask, alternative to specifying xy pairs, may be useful if roipoly is used\n');
fprintf('''yCoord'',[ y1 y2 ...]               y-coordinates of the beamstop mask, alternative to specifying xy pairs, may be useful if roipoly is used\n');
fprintf('''FilenameValidMask'',<path and filename>  Matlab file with the valid pixel indices ind_valid,\n');
fprintf('                                    default is %s\n',filename_valid_mask);
fprintf('''SaveData'',<0-no,1-yes>             0 for displaying the result without updating the mask, default is %d\n',...
    save_data);
fprintf('''FigNo'',<integer>                   number of the figure in which the result is displayed, default is %d\n',...
    fig_no);
fprintf('\n');
fprintf('A valid pixel mask can be created using the macro prep_valid_mask.\n')
fprintf('You will find a valid pixel mask in %s but you may consider to measure a new one.\n',...
    filename_valid_mask);
fprintf('\n');
