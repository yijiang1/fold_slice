% full_array = add_to_3D_projection(small_array,full_array, positions_offset, indices,add_values, add_atomic, use_MEX)
% add one small 3D block into a large 3D array with offset given by
% positions_offset vector and perform this operation only for slices selected by
% indiced vector 
%
%  Inputs: 
%       **full_array - array to which the small_array will be added / written
%       **small_array - array used to be added to large array 
%       **positions_offset - [Nangles x 2] offset from (1,1) coordinate in pixels
%       for each slice , if provide only [1x2] vector, assume the same
%       offset for each slice 
%       **indices - add only to selected sliced of the full_array
%   *optional*
%       **add_values - (default==true) add values instead of rewritting 
%       **add_atomic - (default==true) add values in atomic way, slow but it allows overlapping regions 
%       **use_MEX -    (use_MEX==true) use fast mex code 
%   *returns*
%        ++full_array or none, results were writted !directly! to the input
%        array full_array, there is not need to take any output if MEX
%        function add_to_3D_projection was used 
%
% Compilation from Matlab:
%     mex -R2018a 'CFLAGS="\$CFLAGS -fopenmp"' LDFLAGS="\$LDFLAGS -fopenmp" add_to_3D_projection_mex.cpp
%   Usage from Matlab:
%  
%   full_array = (rand(1000, 1000, 200, 'single'));
%   small_array = (zeros(500, 500, 100, 'single'));
% 
%   positions_offset = (10*rand(100,2));
%   indices = ([1:100]);  % indices are starting from 1 !! 
%   add_values = true; 
%   add_to_3D_projection(small_array,full_array,positions_offset, indices,add_values);



% *-----------------------------------------------------------------------*
% |                                                                       |
% |  Except where otherwise noted, this work is licensed under a          |
% |  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
% |  International (CC BY-NC-SA 4.0) license.                             |
% |                                                                       |
% |  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
% |                                                                       |
% |       Author: CXS group, PSI                                          |
% *-----------------------------------------------------------------------*
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
% 
% 

function full_array = add_to_3D_projection(small_array,full_array, positions_offset, indices,add_values, add_atomic, use_MEX)
    
    if nargin < 7
        use_MEX = true; 
    end
    if nargin < 6 
        add_atomic = true; 
    end
    if nargin < 5 
       add_values = true;  
    end
    if size(positions_offset,1)==1
        positions_offset = repmat(positions_offset, numel(indices), 1);
    end
    if use_MEX && ~isa(full_array, 'gpuArray') && ~verLessThan('matlab', '9.4') && ~islogical(small_array)  % logical arrays not yet implemented
        %% run fast MEX-based code if possible 
        try
            add_to_3D_projection_mex(small_array,full_array, int32(positions_offset), int32(indices),add_values>0,add_atomic>0);
        catch err 
            % recompile the scripts if needed
            if any(strcmp(err.identifier, { 'MATLAB:UndefinedFunction','MATLAB:mex:ErrInvalidMEXFile'}))
                utils.verbose(0, 'Recompilation of MEX functions ... ')
                path = replace(mfilename('fullpath'), mfilename, ''); 
                mex('-R2018a','-O', 'CFLAGS="\$CFLAGS -fopenmp"', '-O','LDFLAGS="\$LDFLAGS -fopenmp"',[path,'private/add_to_3D_projection_mex.cpp'], '-output', [path, 'private/add_to_3D_projection_mex'])
                add_to_3D_projection_mex(small_array,full_array, int32(positions_offset), int32(indices),add_values>0,add_atomic>0);
            else
               rethrow(err) 
            end
        end
       return 
    end

    %% matlab alternative to the MEX file , (much slower)
    positions_offset = round(positions_offset);
    N_f = size(full_array);
    N_s = size(small_array);

    for ii = 1:size(positions_offset,1)
        jj = min(indices(ii),size(full_array,3)); 
        for i = 1:2
            ind_f{i} = max(1, 1+positions_offset(ii,i)):min(N_f(i),positions_offset(ii,i)+N_s(i));
            ind_s{i} = ((ind_f{i}(1)-positions_offset(ii,i))):(ind_f{i}(end)-positions_offset(ii,i));
        end
        if add_values
            full_array(ind_f{:},jj) = full_array(ind_f{:},jj) + small_array(ind_s{:},min(ii, size(small_array,3)));
        else
            full_array(ind_f{:},jj) = small_array(ind_s{:},min(ii, size(small_array,3))); 
        end
    end

end

