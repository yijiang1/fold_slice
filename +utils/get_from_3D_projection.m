% GET_FROM_3D_PROJECTION add one small 3D block into large 3D array 
%
% small_array = get_from_3D_projection(small_array,full_array, positions_offset, indices)
%
%  Inputs: 
%       **full_array - array from which the small_array will loaded
%       **small_array - empty array for storing the data 
%       **positions_offset - [Nangles x 2] offset from (1,1) coordinate in pixels
%       for each slice , if provide only [1x2] vector, assume the same
%       offset for each slice 
%       **indices - add only to selected sliced of the full_array
%       *optional*
%       **use_MEX -    (use_MEX==true) use fast mex code 
%   *returns*
%        ++small_array or none, results were writted !directly! to the input
%               array small_array, there is not need to take any output if MEX
%               function add_to_3D_projection was used 
%
%  Compilation from Matlab:
%        mex -R2018a 'CFLAGS="\$CFLAGS -fopenmp"' LDFLAGS="\$LDFLAGS -fopenmp" get_from_3D_projection_mex.cpp
%   Usage from Matlab:
%  
%   full_array = (randn(1000, 1000, 1, 'single'));
%   small_array = (ones(500, 500, 100, 'single'));
% 
%   positions_offset = int32([1:100; 1:100])';
%   indices = int32([1:100]);  % indices are starting from 1 !! 
%   tic; get_from_3D_projection_mex(small_array,full_array,positions_offset,indices); toc



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
% 



function small_array = get_from_3D_projection(small_array, full_array, positions_offset, indices, use_MEX)

    Np_f = size(full_array);
    Np_s = size(small_array);

    if nargin < 4
        indices = 1:Np_f(3);  
    end
    if nargin < 5
        use_MEX = true; 
    end
    if size(positions_offset,1)==1
        positions_offset = repmat(positions_offset, size(small_array,3), 1);
    end
    positions_offset = int32(positions_offset);
    indices= int32(indices); 
    
    if use_MEX && ~isa(full_array, 'gpuArray') && ~verLessThan('matlab', '9.4') && ~islogical(small_array)  % logical arrays not yet implemented
        %% run fast MEX-based code if possible 
        try
            get_from_3D_projection_mex(small_array,full_array, positions_offset, indices)
        catch err 
            % recompile the scripts if needed
            if any(strcmp(err.identifier, { 'MATLAB:UndefinedFunction','MATLAB:mex:ErrInvalidMEXFile'}))
                utils.verbose(0, 'Recompilation of MEX functions ... ')
                path = replace(mfilename('fullpath'), mfilename, ''); 
                mex('-R2018a','-O', 'CFLAGS="\$CFLAGS -fopenmp"', '-O','LDFLAGS="\$LDFLAGS -fopenmp"',[path,'private/get_from_3D_projection_mex.cpp'], '-output', [path, 'private/get_from_3D_projection_mex'])
                get_from_3D_projection_mex(small_array,full_array, positions_offset, indices)
            else
               rethrow(err) 
            end
        end
       return 
    end


    
    %% simple matlab based version that may be too slow 
    for ii = 1:size(positions_offset,1)
        jj = min(indices(ii),size(full_array,3)); 
        for i = 1:2
            % limit to the region inside full_array
            ind_f{i} = positions_offset(ii,i)+int32(1:Np_s(i));
            ind_f{i} = max(1,1+positions_offset(ii,i)):min(positions_offset(ii,i)+Np_s(i),Np_f(i));
            % adjust size of the small matrix to correspond
            ind_s{i} = ((ind_f{i}(1)-positions_offset(ii,i))):(ind_f{i}(end)-positions_offset(ii,i));
        end
        small_array(ind_s{:},ii) = full_array(ind_f{:},jj) ;
    end   
end
