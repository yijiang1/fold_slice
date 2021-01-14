%UPDATE_PARAM Use update_param to update a struct with another struct and /
% or fields and their corresponding values.
%   
%   Update struct with another struct:
%     p = update_param(p, p_updated);
%
%   Update struct with another struct and fields:
%     p = update_param(p, p_updated, 'windowautopos', 1, 'extrastringtitle', 'Final recon');
%   
%   Update struct with fields:
%       p = update_param(p, 'windowautopos', 1, 'extrastringtitle','Final recon');
%
%   Use field force_update to select different update types:
%       force_update = 0 ...        add fields but don't overwrite them
%       force_update = 1 ...        overwrite fields but leave protected
%                                   fields unchanged (default)
%       force_update = 2 ...        overwrite everything
%
%       p = update_param(p, 'windowautopos', 1, 'force_update', 0);
%

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
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

function [ p ] = update_param(p, varargin)
import utils.param_protect_field

% get protected fields
prot_fn = param_protect_field(); 

update_struct = [];
if ~isempty(varargin)
    % check if first varargin is struct
    if isstruct(varargin{1})
        update_struct = varargin{1};
        % load additional params
        vararg = cell(0,0);
        if length(varargin)>2
            for ind = 2:2:length(varargin)
                vararg{end+1} = varargin{ind};
                vararg{end+1} = varargin{ind+1};
            end
        end
    else
        % if first argument is not a struct, load params
        vararg = cell(0,0);
        if length(varargin)>1
            for ind = 1:2:length(varargin)
                vararg{end+1} = varargin{ind};
                vararg{end+1} = varargin{ind+1};
            end
        end
    end
end
force_update = 1; 
if ~isempty(update_struct)
    if isfield(update_struct, 'force_update')
        force_update = update_struct.force_update;
    end
end
if ~isempty(vararg)
    for ind = 1:2:length(vararg)
        if strcmp(vararg{ind}, 'force_update')
            force_update = vararg{ind+1};
        end
    end
end

% update struct with struct
if ~isempty(update_struct)
    fn_struct = fieldnames(update_struct);
    for ii=1:size(fn_struct,1)
        if ~isfield(p, fn_struct{ii}) || force_update > 0
            if ~isfield(p, (fn_struct{ii}))
                if strcmp(fn_struct{ii}, 'force_update')
                    continue;
                end
            elseif ~isempty(prot_fn)
                if ismember(fn_struct{ii},prot_fn) && force_update < 2
                    error('You are trying to update the protected field ''%s''! If you are sure, use force_update=2.', fn_struct{ii})
                end
            end
            p.(fn_struct{ii}) = update_struct.(fn_struct{ii});
        end
    end
end

% update struct with fieldnames
if ~isempty(vararg)
    for ind = 1:2:length(vararg)
        if ~isfield(p, vararg{ind})
            if strcmp(vararg{ind}, 'force_update')
                continue;
            end
        elseif ~isempty(prot_fn)
            if ismember(vararg{ind},prot_fn) && force_update < 2
                error('You are trying to update the protected field ''%s''! If you are sure, use force_update=2.', vararg{ind})
            end
        end
        p.(vararg{ind}) = vararg{ind+1};
    end
end

end

