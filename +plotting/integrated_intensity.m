% Call function without arguments for instructions on how to use it

% Filename: $RCSfile: integrated_intensity.m,v $
%
% $Revision: 1.1 $  $Date: 2009/04/08 16:46:40 $
% $Author:  $
% $Tag: $
%
% Description:
% sum the intensity of frames
%
% Note:
% Call without arguments for a brief help text.
%
% Dependencies: 
% - image_read
%
% history:
%
% September 4th 2009: 
% use find_files rather than dir to find the files
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

function [int_int] = integrated_intensity(filename_mask, varargin)
import beamline.pilatus_valid_pixel_roi
import io.image_read
import utils.find_files

% initialize return arguments
valid_mask = struct('indices',[], 'framesize',[]);

% set default values for the variable input arguments:
% filename for loading and saving the valid pixel mask
filename_valid_mask = '~/Data10/analysis/data/pilatus_valid_mask.mat';
% display result in this figure
fig_no = 2;

% check minimum number of input arguments
if (nargin < 1)
    fprintf('\nUsage:\n');
    fprintf('[int_int]=%s(filename_mask [[,<name>,<value>]...]);\n',mfilename);
    fprintf('''FigNo'',<integer>                   number of the figure in which the result is displayed, default is %d\n',...
        fig_no);
    fprintf('\n');
    fprintf('Examples:\n');
    fprintf('[int_int]=%s(''~/Data10/pilatus/air_scattering/*.cbf'');\n',...
        mfilename);

    error('At least the filename mask has to be specified as input parameter.');
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
            no_of_in_arg = no_of_in_arg -1 + length(varargin);
        end
    end
end

% check number of input arguments
if (rem(no_of_in_arg,2) ~= 1)
    error('The optional parameters have to be specified as ''name'',''value'' pairs');
end


% parse the variable input arguments:
% initialize the list of unhandled parameters
vararg_remain = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'FigNo' 
            fig_no = value;
        otherwise
            vararg_remain{end+1} = name; %#ok<AGROW>
            vararg_remain{end+1} = value; %#ok<AGROW>
    end
end

vararg_remain{end+1} = 'UnhandledParError';
vararg_remain{end+1} = 0;


% set some default values for the plot window
set(0, 'DefaultAxesfontsize', 12);
set(0, 'DefaultAxeslinewidth', 1, 'DefaultAxesfontsize', 12);
set(0, 'DefaultLinelinewidth', 1);

% get all matching filenames
[data_dir,fnames,vararg_remain] = ...
    find_files( filename_mask, vararg_remain );

if (length(fnames) < 1)
    error('No matching files found for %s.\n',search_mask);
end

fprintf('loading the valid pixel mask %s\n',filename_valid_mask);
load(filename_valid_mask);


% process the frames
int_int = zeros(1,length(fnames));
fprintf('data directory is %s\n',data_dir);
for (f_ind=1:length(fnames)) 
    fprintf('%3d/%3d: reading %s%s\n',f_ind,length(fnames),...
        data_dir,fnames(f_ind).name);
    [frame] = image_read([data_dir fnames(f_ind).name ],vararg_remain);
    if (f_ind == 1)
        framesize1 = size(frame.data,1);
        framesize2 = size(frame.data,2);
        framesize = framesize1 * framesize2;
        
        ind_valid = 1:numel(frame.data);
        if (isstruct(valid_mask))
            % cut out the region of interest from the valid pixel mask
            valid_mask_cut = pilatus_valid_pixel_roi(valid_mask,...
                'RoiSize',frame.img_full_size{1});
            % from this cut out the region of interest, 
            % if specified
            valid_mask_cut = pilatus_valid_pixel_roi(valid_mask_cut,vararg_remain);

            ind_valid = valid_mask_cut.indices;
        end
    end

    % check that the file have identical dimensions
    if ((framesize1 ~= size(frame.data,1)) || ...
        (framesize2 ~= size(frame.data,2)))
            error('The previous file(s) had %d x %d pixels, this frame has %d x %d pixels',...
            framesize1,framesize2,size(frame.data,1),size(frame.data,2));
    end
    
    % sum the valid pixels
    int_int(f_ind) = sum(frame.data(ind_valid));
    
end


% plot the result
if (fig_no > 0)
    figure(fig_no);
    plot(int_int);
    title_str = [ strrep(filename_mask,'_','\_') ': average int. int. = ' ...
        num2str(mean(int_int),'%.4e') ];
    title(title_str);
end
