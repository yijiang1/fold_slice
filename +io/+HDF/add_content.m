%ADD_CONTENT write matlab structure to H5 file

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

function add_content(data, gid, plist, comp, overwrite)
import io.HDF.*

fn = fieldnames(data);
for jj=1:length(fn)
    if isstruct(data.(fn{jj}))
        fn_sub = fieldnames(data.(fn{jj}));
        if length(fn_sub) <= 2 && isfield(data.(fn{jj}), 'Value')
            % found dataset
            if isfield(data.(fn{jj}), 'Attributes')
                Attributes.MATLAB_class = class(data.(fn{jj}).Value);
                write_dataset(data.(fn{jj}).Value, gid, fn{jj}, plist, comp, overwrite, data.(fn{jj}).Attributes);
            else
                Attributes.MATLAB_class = class(data.(fn{jj}).Value);
                write_dataset(data.(fn{jj}).Value, gid, fn{jj}, plist, comp, overwrite, Attributes);
            end
        elseif strcmpi(fn{jj}, 'Attributes')
            % add attributes to group
            fn_attr = fieldnames(data.Attributes);
            for ii=1:length(fn_attr)
                write_attribute(gid, data.Attributes.(fn_attr{ii}), fn_attr{ii});
            end
        else
            % found group
            gid_new = add_groups(gid, fn{jj}, plist, true);
            if length(data.(fn{jj})) == 1
                add_content(data.(fn{jj}), gid_new, plist, comp, overwrite)
            else
                fn_names = cell(1,length(data.(fn{jj})));
                for ii=1:length(data.(fn{jj}))
                    fn_names{ii} = sprintf([fn{jj} '_%d'],ii-1);
                    gid_new_sub = add_groups(gid_new, fn_names{ii}, plist, true);
                    add_content(data.(fn{jj})(ii), gid_new_sub, plist, comp, overwrite);
                end
                write_attribute(gid_new, 'structure array', 'MATLAB_class');

            end
        end
    else
        % append dataset
        Attributes.MATLAB_class = class(data.(fn{jj}));
        write_dataset(data.(fn{jj}), gid, fn{jj}, plist, comp, overwrite, Attributes);
    end
end
        
end
