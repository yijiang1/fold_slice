% Call function without arguments for instructions on how to use it

% Filename: $RCSfile: image_orient.m,v $
%
% $Revision: 1.5 $  $Date: 2014/04/11 12:47:09 $
% $Author:  $
% $Tag: $
%
% Description:
% Macro for mirroring or rotating images or stacks of images. 
%
% Note:
% Call without arguments for a brief help text.
% In case of image stacks the first two dimensions are treated as the image
% dimensions. 
%
% Dependencies: 
% ---
%
% history:
%
% January 16th 2009:
% return the complete structure rather than just the data array
%
% June 19th 2008:
% change call to complete frame rather than data only
%
% May 27th 2008: 
% 1st version based on orientm.m by Tilman Donath

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

function [frame_out,vararg_remain] = image_orient(frame, varargin)
import io.*
import utils.default_parameter_value

% check minimum number of input arguments
if (nargin < 1)
    image_orient_help(mfilename);
    error('At least one input parameter has to be specified.');
end

if (ndims(frame.data) < 2)
    error('The input data array must have at least two dimensions.');
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
    error('The optional parameters have to be specified as ''name'',''value'' pairs');
end


% set default values for the variable input arguments:
orient_extension = frame.extension{1};
do_transpose = default_parameter_value(mfilename,'Transpose');
do_fliplr = default_parameter_value(mfilename,'FlipLR');
do_flipud = default_parameter_value(mfilename,'FlipUD');
invert_orientation = 0;

% parse the variable input arguments
vararg_remain = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'Transpose' 
            do_transpose = value;
        case 'FlipLR' 
            do_fliplr = value;
        case 'FlipUD' 
            do_flipud = value;
        case 'Orientation' 
            if (length(value) ~= 3)
                error('Invalid Orientation parameter of length %d',length(value));
            end
            do_transpose = value(1);
            do_fliplr = value(2);
            do_flipud = value(3);
        case 'OrientExtension' 
            orient_extension = value;
        case 'OrientByExtension' 
            if (value)
                orient_vec = image_default_orientation(frame.header,orient_extension);
                do_transpose = orient_vec(1);
                do_fliplr = orient_vec(2);
                do_flipud = orient_vec(3);
            end
        case 'InvertOrientation'
            invert_orientation = value;
        otherwise
            vararg_remain{end+1} = name;
            vararg_remain{end+1} = value;
    end
end


% initialize return arguments
frame_out = frame;

if (invert_orientation)
    % apply orientation modifications in inverse order, e.g., to revert the
    % original orientation prior to writing a file
    if (do_flipud)
        frame_out.data = flip(frame_out.data,1);
    end
    
    if (do_fliplr)
        frame_out.data = flip(frame_out.data,2);
    end
end

if (do_transpose)
    dim_order = 1:ndims(frame_out.data);
    dim_order(1) = 2;
    dim_order(2) = 1;
    frame_out.data = permute(frame_out.data,dim_order);
end

if (~invert_orientation)
    if (do_fliplr)
        frame_out.data = flip(frame_out.data,2);
    end

    if (do_flipud)
        frame_out.data = flip(frame_out.data,1);
    end
end
