%HDF5_CP_FILE copy HDF files
%   orig_filename...        source file
%   duplicate_filename...   target file
%
%   *optional*              given as name/value pair
%   groups...               groups to copy; either string or cell of
%                           strings; default: everything in root 
%   copy_type...            'deep', 'normal' or 'shallow' copy; 
%                           'shallow' creates external links in target file; 
%                           'normal' is similar to linux 'cp' command; 
%                           'deep' dereferences all internal and external links; 
%                           default: 'shallow'
%
%   EXAMPLES:
%       hdf5_cp_file('./test.h5', './test_new.h5')
%       hdf5_cp_file('./test.h5', './test_new.h5', 'copy_type', 'deep');
%   
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

function hdf5_cp_file(orig_filename, duplicate_filename, varargin)
import io.HDF.*
% take care of input arguments
groups = [];
copy_type = 'shallow';

% parse the variable input arguments vararg = cell(0,0);
if ~isempty(varargin)
    for ind = 1:2:length(varargin)
        name = varargin{ind};
        value = varargin{ind+1};
        switch lower(name)
            case 'groups'
                groups = value;
            case 'copy_type'
                copy_type = value;
        end
    end
end

switch copy_type
    case 'shallow'
        if isempty(groups)
            % if no groups are specified, use h5info to get all datasets and groups
            % from root
            h = h5info(orig_filename, '/');
            lng = length(h.Groups);
            lnd = length(h.Datasets);
            lna = length(h.Attributes);
            
            groups = cell([1 lng+lnd]);
            attributes = [];
            
            for ii=1:lng
                groups{ii} = h.Groups(ii).Name;
            end
            for ii=1:lnd
                groups{ii+lng} = h.Datasets(ii).Name;
            end
            for ii=1:lna
                attributes.(h.Attributes(ii).Name) = h.Attributes(ii).Value;
                if iscell(h.Attributes(ii).Value)
                    attributes.(h.Attributes(ii).Name) = attributes.(h.Attributes(ii).Name){1};
                end
            end
        else
            attributes = [];
        end
        
        
        s = [];
        if iscell(groups)
            for ii=1:length(groups)
                subgrps = strsplit(rm_delimiter(groups{ii}), '/');
                s = setfield(s, subgrps{:}, ['ext:' orig_filename ':' groups{ii}]);
            end
        else
            s.groups = ['ext:' orig_filename ':' groups];
        end
        
        % append attributes
        if ~isempty(attributes)
            s.Attributes = attributes;
        end
        
        save2hdf5(duplicate_filename, s, 'overwrite', true, 'iscopy', true);
        
    case 'deep'
        s = io.HDF.hdf5_load(orig_filename, '-ca');
        save2hdf5(duplicate_filename, s, 'overwrite', true, 'iscopy', true);
        
    case 'normal'
        copyfile(orig_filename, duplicate_filename)
        
    otherwise
        error('Unknown copy type!')
end

end
