% Call function without arguments for instructions on how to use it

% Filename: $RCSfile: pilatus_valid_pixel_roi.m,v $
%
% $Revision: 1.5 $  $Date: 2011/05/19 16:44:24 $
% $Author:  $
% $Tag: $
%
% Description:
% cut out of the valid pixel mask for the full detector the one for the
% current region of interest
%
% Note:
% The location of the ROI is not stored in the data files. It is deduced
% from the known readoiut modes. So far only the 1x2 module mode is
% implemented. 
%
% Dependencies: 
% none
%
% history:
%
% January 31st 2009: add arbitrary ROIs via named parameters
%
% May 9th 2008: 1st version

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

function [valid_mask] = pilatus_valid_pixel_roi(valid_mask,varargin)

% sub-detector readout size
roi_size = [];

% alternatively:
% from/to row 0 means all rows
row_from = 0;
row_to = 0;
% from/to column 0 means all lines
column_from = 0;
column_to = 0;

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
if ((nargin < 1) || (rem(no_of_in_arg,2) ~= 1))
    fprintf('Usage:\n')
    fprintf('%s(valid_mask,[[<name>,<value>], ...]);\n',mfilename);
    fprintf('The name value pairs are:\n');
    fprintf('''RoiSize'',[<size-y> <size-x>]        size of sub-detector readout ROIs\n');
    fprintf('Alternatively:\n');
    fprintf('''RowFrom'',<0-max>                    region of interest definition, 0 for full frame\n');
    fprintf('''RowTo'',<0-max>                      region of interest definition, 0 for full frame\n');
    fprintf('''ColumnFrom'',<0-max>                 region of interest definition, 0 for full frame\n');
    fprintf('''ColumnTo'',<0-max>                   region of interest definition, 0 for full frame\n');
    fprintf('''ROI'',s[ <ColumnFrom> <RowFrom> <ColumnTo> <RowTo> ]\n');
    fprintf('                                     region of interest definition of all four coordinates together\n');
end


% parse the variable input arguments
vararg = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'ROI'
            if (length(value) ~= 4)
                error('The ROI parameter needs a vector of length four as argument.');
            end
            column_from = value(1);
            row_from = value(2);
            column_to = value(3);
            row_to = value(4);
        case 'RowFrom' 
            row_from = value;
            roi_size = [];
        case 'RowTo' 
            row_to = value;
            roi_size = [];
        case 'ColumnFrom' 
            column_from = value;
            roi_size = [];
        case 'ColumnTo' 
            column_to = value;
            roi_size = [];
        case 'RoiSize' 
            roi_size = value;
        otherwise
            vararg{end+1} = name; %#ok<AGROW>
            vararg{end+1} = value; %#ok<AGROW>
    end
end


% if (valid_mask.framesize ~= [1679 1475])
%     error('can not handle the source size (%d,%d)',...
%         valid_mask.framesize(1),valid_mask.framesize(2));
% end


if (isempty(roi_size))
    if (row_from < 1)
        row_from = 1;
    end
    if (row_to < 1)
        row_to = valid_mask.framesize(1);
    end
    if (column_from < 1)
        column_from = 1;
    end
    if (column_to < 1)
        column_to = valid_mask.framesize(2);
    end
    % nothing to do
    if ((row_from == 1) && (column_from == 1) && ...
        (row_to == valid_mask.framesize(1)) && (column_to == valid_mask.framesize(2)))
        return;
    end
    
    x_from = column_from;
    y_from = row_from;
    roi_size = [ row_to-row_from+1 column_to-column_from+1 ];
else
    % nothing to do
    if (valid_mask.framesize == roi_size)
        return;
    end

    % determine the location of the ROI from its size via the known modi
    x_from = 0;
    y_from = 0;
    % two modules 
    if (roi_size == [407 487])
        x_from =  495;
        y_from =  637;
    end
    if (roi_size == [831 1475])
        x_from =  1;
        y_from =  425;
    end
    if (roi_size == [831 981])
        x_from = 495;
        y_from = 425;
    end
    if ((x_from == 0) || (y_from == 0))
        error('can not handle the ROI size (%d,%d)',roi_size(1),roi_size(2));
    end
end




% create the full valid pixel mask
vpm = zeros(valid_mask.framesize);
vpm(valid_mask.indices) = 1;
% cut out the region of interest
vpm = vpm(y_from:(y_from+roi_size(1)-1), ...
          x_from:(x_from+roi_size(2)-1));

% return the indices of valid pixels within this ROI
valid_mask.indices = find(vpm == 1);
valid_mask.framesize = size(vpm);
