% [outstr] = char_to_cellstr(inchars,nl_only)
% Convert an array of text to a cell array of lines.

% Filename: $RCSfile: char_to_cellstr.m,v $
%
% $Revision: 1.4 $  $Date: 2014/04/11 10:57:20 $
% $Author:  $
% $Tag: $
%
% Description:
% Convert an array of text to a cell array of lines.
%
% Note:
% Used for making file headers accessible. 
%
% Dependencies:
% none
%
%
% history:
%
% June 22nd 2008: bug fix for fliread adding the nl_only 
% parameter, to be replaced by named parameter later on
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

function [outstr] = char_to_cellstr(inchars,nl_only)

if (nargin < 2)
    nl_only = 0;
end

% get positions of end-of-line signatures
eol_ind = regexp(inchars,'\r\n');
eol_offs = 1;
if ((length(eol_ind) < 1) || (nl_only))
    eol_ind = regexp(inchars,'\n');
    eol_offs = 0;
end
if (length(eol_ind) < 1)
    eol_ind = length(inchars) +1;
end
if (length(eol_ind) < 1)
    outstr = [];
    return;
end

% dimension return array with number of lines
outstr = cell(length(eol_ind),1);

% copy the lines to the return array, not suppressing empty lines
start_pos = 1;
ind_out = 1;
for (ind = 1:length(eol_ind))
    end_pos = eol_ind(ind) -1;
    % cut off trailing spaces
    while ((end_pos >= start_pos) && (inchars(end_pos) == ' '))
        end_pos = end_pos -1;
    end
    % store non-empty strings
    if (end_pos >= start_pos)
        outstr{ind_out} = inchars(start_pos:end_pos);
    else
        outstr{ind_out} = '';
    end
    ind_out = ind_out +1;
    
    start_pos = eol_ind(ind) +1 + eol_offs;
    ind = ind +1;
end

% resize cell array in case of empty lines
if (ind_out <= length(eol_ind))
    outstr = outstr(1:(ind_out-1));
end

