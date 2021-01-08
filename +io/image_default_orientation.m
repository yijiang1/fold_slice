% [orient_vec] = image_default_orientation(header, extension, varargin)
% Determine default orientation for an image_orient.m call based on the
% file extension

% Filename: $RCSfile: image_default_orientation.m,v $
%
% $Revision: 1.11 $  $Date: 2013/03/23 15:01:09 $
% $Author:  $
% $Tag: $
%
% Description:
% Determine default orientation for an image_orient.m call based on the
% file extension
%
% Note:
% Call without arguments for a brief help text.
%
% Dependencies: 
% ---
%
% history:
%
% September 30th 2010:
% add orientation for HDF5 files
%
% November 19th 2008:
% add .mat files
%
% August 28th 2008:
% add orientation for extension .dat
%
% July 17th 2008:
% add mar extension, raw default orientation changed before
%
% June 19th 2008:
% add header to call parameters
%
% June 10th 2008: 
% 1st version

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

function [orient_vec] = ...
   image_default_orientation(header, extension, varargin)
import io.*

% check minimum number of input arguments
if (nargin < 1)
    fprintf('Usage:\n')
    fprintf('[orientation_vector]=%s(extension);\n',...
        m_file_name);
    fprintf('The vector contains three values which can be 0 or 1 for transpose, flip-left-right, flip-up-down\n');
    error('At least one input parameter has to be specified.');
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
% vararg_remain = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
%     value = varargin{ind+1};
    switch name
        otherwise
            error('Do not know how to handle parameter %s\n',name);
%             vararg_remain{end+1} = name;
%             vararg_remain{end+1} = value;
    end
end


% image_read calls this function prior to converting the single frame
% to a cell array of frames
if (~iscell(header))
    fprintf('Warning (%s): header is not a cell array\n',mfilename);
    fprintf('If this is an image_spec call then please report to Oliver:\n')
    whos header
    header
else
    if ((~isempty(header)) && (iscell(header{1})))
        header = header{1};
    end
end


% set the default orientation as a function of the filename extension
switch extension
    case 'dat'
        orient_vec = [ 0 0 0 ];
    case 'edf'
        orient_vec = [ 1 1 1 ];
    case 'cbf'
        orient_vec = [ 1 1 1 ];
    case {'h5', 'hdf5', 'nxs', 'cxs'}
        orient_vec = [ 0 0 1 ];        
    case {'tif', 'tiff'}
        if (strcmp(header{1}(1:5),'Andor'))
            orient_vec = [ 1 0 0 ];
        else
            orient_vec = [ 0 0 0 ];
        end
    case {'mar','mccd'}
        orient_vec = [ 1 0 1 ];
    case 'mat'
        orient_vec = [ 0 0 0 ];
    case 'raw'
        % FLI CCD at ICON
%        orient_vec = [ 0 1 0 ];
        % FLI CCD at laser setup
        orient_vec = [ 1 0 1 ];
    case 'spe'
        orient_vec = [ 0 0 0 ];
    otherwise
        error([ 'unknown extension ''' extension '''' ]);
end
