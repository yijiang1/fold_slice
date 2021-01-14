% Function verbose(level, message)
%
% verbose : shows current verbose level
% verbose(n): set the level to n
% verbose(n, message): displays the message if n <= current verbose level
% 
% verbose(n, message, v1, v2, ...) behaves like sprintf(message, v1, v2, ...)
%
% Suggested levels:
%  1: important information
%  2: general information
%  3: debugging

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

function varargout = verbose(varargin)

persistent verbose_level
persistent verbose_prefix

% The output function to use (this could be configurable)
do_output = @output_print_fct;

if isempty(verbose_level)
    verbose_level = 1;
end

if nargin == 0
    varargout{1} = verbose_level;
    return
end

if nargin == 1
    if ischar(varargin{1})
        verbose_level = str2num(varargin{1});
    elseif isstruct(varargin{1})
        verbose_prefix = varargin{1}.prefix;
    else
        verbose_level = varargin{1};
    end
    return
end

if nargin >= 2
    if verbose_level >= varargin{1}
        do_output(varargin{1}, sprintf(varargin{2},varargin{3:end}), verbose_level, verbose_prefix);
    end
end    

% Two examples of output functions.
function output_print(level, str)
disp(str);

function output_print_fct(level, str, verbose_level, verbose_prefix)
if verbose_level > 3
    [st, ~] = dbstack(2);
    if ~isempty(st) 
        str = sprintf('%s [%d] : %s', st(1).name, st(1).line, str);
    else
        str = sprintf('[root] : %s', str);
    end
elseif ~isempty(verbose_prefix)
    str = sprintf('[%s] : %s', verbose_prefix, str);
end
disp(str);


