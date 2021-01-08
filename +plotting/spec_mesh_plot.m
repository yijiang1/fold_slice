% function h = plotting.spec_mesh_plot(spec_dat_file, scan, ['counter',counter_name])
%
%SPEC_MESH_PLOT Plot the counter values from a mesh scan. In the future to
%   be generalized to other types of scans
%
% inputs
%
% ** spec_dat_file  Filename and path of the SPEC dat file. A base path can
%                   also be given and the code will search for the file
% ** scan           Scan number
%
%*optional*
% ** 'counter',counter_name     counter_name is the name of the SPEC
%                               counter, by default it will be 'bpm4i'
%
% returns
% ++ h              Figure handle
%
%
% see also: 
%
% EXAMPLES:
%   plotting.spec_mesh_plot('~/Data10', 10)

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2019 by Paul Scherrer Institute (http://www.psi.ch)    |
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


function h = spec_mesh_plot(spec_dat_file, scan, varargin)


% set default values for the variable input arguments:
counter = 'bpm4i';

% check minimum number of input arguments
if (nargin < 2)
    error('At least the spec dat filename and scan number have to be specified as input parameters.');
end

% accept cell array with name/value pairs as well
no_of_in_arg = nargin;

% check number of input arguments
if (rem(no_of_in_arg,2) ~= 0)
    error('The optional parameters have to be specified as ''name'',''value'' pairs');
end


% parse the variable input arguments
vararg_remain = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch lower(name)
        case 'counter' 
            counter = value;
        otherwise
            vararg_remain{end+1} = name; %#ok<AGROW>
            vararg_remain{end+1} = value; %#ok<AGROW>
    end
end

S = io.spec_read(spec_dat_file,'ScanNr',scan);
scanstring = strsplit(S.S);
scantype = scanstring{3};

switch lower(scantype)
    case 'mesh'
        name_axis_fast = scanstring{4};
        N_fast = str2num(scanstring{7})+1;
        name_axis_slow = scanstring{8};
        N_fast = str2num(scanstring{7})+1;
        N_slow = str2num(scanstring{11})+1;
    otherwise
        error(['Scan type ' lower(scantype) ' is not defined for this function'])
end
    

%%%%  HASTA AQUI VOY %%%%

data      = reshape(S.(counter),N_fast,N_slow).';
axis_fast = reshape(S.(name_axis_fast),N_fast,N_slow).';
axis_slow = reshape(S.(name_axis_slow),N_fast,N_slow).';
if nargout > 0
    h = figure(1);
else
    figure(1)
end
imagesc(axis_fast(1,:),axis_slow(:,1),data)
title(S.S);
xlabel(name_axis_fast)
ylabel(name_axis_slow)
axis xy

end

