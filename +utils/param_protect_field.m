%%PARAM_PROTECT_FIELD
%   param_protect_field(param)... check for protected field; returns
%   boolean
%
%   accepts struct or string as input
%
%   param_protect_field()...    return protected fields
%
%   param_protect_field(param, 'p')...  add param to protected fields
%
%   param_protect_field(param, 'r')...  remove param from protected fields
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

function [varout] = param_protect_field(varargin)

persistent prot_field

if nargin == 0
    if isempty(prot_field)
        varout = [];
    elseif isempty(fieldnames(prot_field))
        varout = [];
    else
        varout = fieldnames(prot_field);
    end
    return
end

if nargin == 1 
    
    if isstruct(varargin{1})
        fn = fieldnames(varargin{1});
        for ii=1:length(fn)
            varout{ii} = isfield(prot_field, fn{ii});
        end
    elseif ischar(varargin{1})
        varout{1} = isfield(prot_field,varargin{1});
    end
    
elseif nargin == 2
    
    if strcmp(varargin{2},'p')
        % protect fields
        prot_field.(varargin{1}) = true;
        
    elseif strcmp(varargin{2}, 'r')
        % remove protected fields
        if isstruct(varargin{1})
            fn = fieldnames(varargin{1});
            for ii=1:length(fn)
                try
                    prot_field = rmfield(prot_field,fn{ii});
                catch
                    fprintf('Could not find protected field %s\n', fn{ii});
                end
            end
            
        elseif ischar(varargin{1})            
            try
                prot_field = rmfield(prot_field,varargin{1});
            catch
                fprintf('Could not find protected field %s\n', varargin{1});
            end
        end
    else
        error('Unknown second argument %s. Please use ''p'' to protect and ''r'' to remove %s from protected fields.', varargin{2}, varargin{1})
    end
    
end

if isempty(prot_field)
    varout = [];
end


end

