% [] = image_orient_help(m_file_name,varargin)
% parameter help for image_orient and calling functions

% Filename: $RCSfile: image_orient_help.m,v $
%
% $Revision: 1.3 $  $Date: 2009/01/16 15:30:29 $
% $Author:  $
% $Tag: $
%
% Description:
% parameter help for image_orient and calling functions
%
% Note:
% none
%
% Dependencies:
% none
%
%
% history:
%
% May 27th 2008: 1st version

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

function [] = image_orient_help(m_file_name,varargin)
import io.*

% check minimum number of input arguments
if (nargin < 1)
    error('At least the m-file name has to be specified as input parameter.');
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


% parse the variable input arguments
parameters_only = 0;
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'ParametersOnly' 
            parameters_only = value;
        otherwise
            error('unknown parameter %s',name);
    end
end

if (~parameters_only)
    fprintf('Usage:\n')
    fprintf('[data_out]=image_orient(<frame> [[,<name>,<value>] ...]);\n');
    fprintf('The optional <name>,<value> pairs are:\n');
end
fprintf('''OrientExtension'',<''extension''>      file name extension to determine the orientation from\n');
fprintf('''OrientByExtension'',<0-no,1-yes>     use the default orientation for this file type, default yes,\n');
fprintf('                                     superseeded by following orientation parameters (i.e., parameter order matters)\n');
fprintf('''Transpose'',<0-no,1-yes>             mirror at the diagonal\n');
fprintf('''FlipUD'',<0-no,1-yes>                mirror about the horizontal axis\n');
fprintf('''FlipLR'',<0-no,1-yes>                mirror about the vertical axis\n');
fprintf('''Orientation'',<[<Transpose> <FlipLR> <FlipUD>]>\n');
fprintf('                                     specify for each of the three parameters 0 (no) or 1 (yes)\n');
if (~parameters_only)
    fprintf('\n');
    fprintf('Examples:\n');
    fprintf('[data_out]=%s(data, ''Orientation'',[1 1 0]);\n',mfilename);
    fprintf('[data_out]=%s(data, ''Transpose'',1, ''FlipLR'',1);\n',mfilename);
    fprintf('These two calls are equivalent.\n');
end
