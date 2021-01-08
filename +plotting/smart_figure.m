%SMART_FIGURE update figure without stealing focus
% fig = plotting.smart_figure(varargin)
%
% **varargin        same inputs as for MATLAB's figure function
% 
% returns:
% ++ fig            figure handle
% 
%   EXAMPLE:
%       fig = plotting.smart_figure(1);
%
% see also: plotting.imagesc3D 

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


function [varargout] = smart_figure(varargin)

if nargin ~= 1
    % use MATLAB's figure function for name-value pairs
    fig = figure(varargin{:});
else
    try
        % check if figure already exists and retrieve ID
        id = varargin{1};
        figList = get(groot, 'Children');
        figID = find([figList==id]);
        
        % if the figure exists, set it as the current figure
        if ~isempty(figID)
            fig = figList(figID);
            set(groot, 'CurrentFigure', fig)
        else
            % figure does not exist - create new one
            fig = figure(id);
        end
    catch
        % ... just in case
        fig = figure(varargin{:});
    end
end

% return figure handle if needed
if nargout>0
    varargout{1} = fig;
end

end


