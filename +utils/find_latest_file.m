%FIND_LATEST_FILE find latest file in given directory and return path
%   
%   *optional input*
%   path...                     search path; default './'
%   file mask...                limit results to a specific name or file
%                               extension; default none
%   offset...                   take latest-offset; default 0
%
%
%   EXAMPLE:
%       out = find_latest_file;
%       out = find_latest_file('../analysis');
%       out = find_latest_file('../analysis', '*.h5');
%       out = find_latest_file('../analysis', '*recons*.h5');
%       out = find_latest_file('../analysis', {*recons*.h5, *recons*.mat});
%       out = find_latest_file('../analysis', '*.h5', -2);
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

function [ out ] = find_latest_file( varargin )

% input
vars = [];

if nargin>0
    vars.path = varargin{1};
else
    vars.path = './';
    vars.name = [];
    vars.offset = 0;
end

if nargin>1
    vars.name = varargin{2};
    vars.offset = 0;
end

if nargin>2
    vars.offset = varargin{3};
end


% make sure that directory exists
if ~isdir(vars.path)
    error('Could not find directory %s', vars.path)
end


% add -name flag if needed
if isempty(vars.name)
    sys_call = sprintf('find %s -type f', vars.path);
else
    if iscell(vars.name)
        nme = ['''' strjoin(vars.name(:), ''' -o -name ''')];
        sys_call = sprintf('find %s -type f \\( -name %s'' \\)', vars.path, nme);
    else
        sys_call = sprintf('find %s -type f -name ''%s''', vars.path, vars.name);
    end
end


% okay, let's go
[~, out] = system([sys_call ' -printf ''%T@ %p\n'' | sort -n | tail ' num2str((abs(vars.offset)*(-1)-1)) '| cut -f2- -d" " | sed -n ''1p''']);

out = out(1:end-1);


end

