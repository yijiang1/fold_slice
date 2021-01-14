% [mask_coord,bmask_ind] = choose_beamstop_mask(filename,varargin)

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

function [mask_coord,bmask_ind] = choose_beamstop_mask(filename,varargin)
import beamline.beamstop_mask
import plotting.image_show

% set default values for the variable input arguments:
% File with mask coordinates
filename_coord = '~/Data10/analysis/data/mask_coordinates.mat';
% border size
border = 3;
% save coordinates
save_coord = 1;
% select corrdinates 
select_points = 1;
% do not read prevously saved coordinates
read_coord = 1;
% start with an empty set of coordinates
mask_coord = [];
bmask_ind = [];
% figure number for display
fig_no_sel = 555 ;
% run 'beamstop_mask' at the end
create_mask = 1;
% Arguments to be passed to imageshow
imageshow_args = {};


% check minimum number of input arguments
if (nargin < 1)
    fprintf('\n');
    fprintf('Usage:\n');
    fprintf('%s(filename,[[<name>,<value>],...]);\n',mfilename);
    fprintf('\n');
    fprintf('The optional <name>,<value> pairs are:\n');
    fprintf('''Border'',                            size of the ''sticky'' edge border, default is %d, 0 to disable\n',border);
    fprintf('''SaveCoord'',<0-no,1-yes>             1 for saving the coordinates, default is %d\n',save_coord);
    fprintf('''FilenameCoord'',<path and filename>  Matlab file with the mask coordinates,\n');
    fprintf('                                     default is %s\n',filename_coord);
    fprintf('''SelectPoints'',<0-no,1-yes>          1 for selecting the points in the image, default is %d\n',select_points);
    fprintf('                                     if 0, points should be either read from the file or \n');
    fprintf('                                     supplied as options for ''beamstop_mask'' function\n');
    fprintf('''ReadCoord'',<0-no,1-yes>             1 for reading the coordinates from the file, default is %d\n',read_coord);
    fprintf('''CreateMask'',<0-no,1-yes>            1 for running ''beamstop_mask'', default is %d\n',create_mask);
    fprintf('''FigNoSel'',<integer>                 number of the figure in which the coordinates are selected, default is %d\n',...
        fig_no_sel);
    fprintf('''ImageShowArgs'', cell                additional parameters to be passed to image_show, default is an empty cell {} \n');
    fprintf('\n');
    fprintf('Additional <name>,<value> pairs recognized by ''beamstop_mask'' can be specified.\n');
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
    display_help(filename_coord,save_coord,create_mask,select_points,border,read_coord,fig_no_sel);
    error('The optional parameters have to be specified as ''name'',''value'' pairs');
end


% parse the variable input arguments
vararg_remain = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'SaveCoord'
            save_coord = value;
        case 'Border'
            border = value;
        case 'SelectPoints'
            select_points = value;
        case 'ReadCoord'
            read_coord = value;
        case 'FilenameCoord'
            filename_coord = value;
        case 'FigNoSel'
            fig_no_sel = value;
        case 'xyCoord'
            mask_coord = value;
        case 'CreateMask'
            create_mask = value;
        case 'ImageShowArgs'
            imageshow_args = value;
        otherwise
            vararg_remain{end+1} = name;
            vararg_remain{end+1} = value;
    end
end





% Read the coordinates from the file
if (read_coord == 1)
    if (exist(filename_coord,'file'))
        load(filename_coord);
        % don't use if there are less than two points in the mask
        % makes it impossible to add new points
        if size(mask_coord,1) < 2
            mask_coord = [];
        end
    else
        fprintf('Mask coordinates file %s was not found.\n',filename_coord);
        fprintf('Continuing with no starting mask.\n');
    end
end


% select/change the mask coordinates in the image
if (select_points == 1)
           
    [qq] = image_show(filename,'FigNo',fig_no_sel, imageshow_args{:});
       
%     h=impoly(gca,mask_coord);
%     mask_coord=getPosition(h);
%     addNewPositionCallback(h,@(pos)eval('mask_coord=pos;'));
%     % wait for the changes while the image is open
%     waitfor(fig_no_sel)
    msgbox({'Instructions:', '1) Create a closed polygon around the beamstop (don''t double-click when you finish)'...
        , '2) Adjust the corners of polygon if needed',...
        '3) Double-click on polygon to finish'},'Choose beamstop mask');
    h = impoly(gca,mask_coord);
    mask_coord = wait(h);
      
    % move points to the edge
    if border 
        im_dim(1) = size(qq.data,2);
        im_dim(2) = size(qq.data,1);

        for jj = 1:size(mask_coord,1)
            for kk = 1:2

                if mask_coord(jj,kk) < border
                    mask_coord(jj,kk) = 0;
                end

                if abs(mask_coord(jj,kk) - im_dim(kk)) < border
                    mask_coord(jj,kk) = im_dim(kk);
                end
            end
        end
    end
    
    % mask coordinates should be integers
    mask_coord = round(mask_coord);
end

% save the mask coordinates file, if specified
if save_coord
    % create a backup of the mask coordinates file
    if (exist(filename_coord,'file'))
        filename_coord_backup = [ filename_coord '.bak' ];
        fprintf('Copying the current mask coordinates file %s to %s\n',filename_coord,...
            filename_coord_backup);
        copyfile(filename_coord,filename_coord_backup);
    end
    fprintf('saving updated mask coordinates file %s\n',filename_coord);
    save(filename_coord,'mask_coord');
end

% run the beamstop_mask function, if specified
if create_mask
    vararg_remain{end+1} = 'xyCoord';
    vararg_remain{end+1} = mask_coord;
    bmask_ind = beamstop_mask(filename,vararg_remain{:},imageshow_args{:});
end







