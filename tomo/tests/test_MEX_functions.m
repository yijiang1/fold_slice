%% simple script to test performace of the MEX accelerated functions 
% utils.get_from_3D_projection and utils.add_to_3D_projection
% the MEX functions use OpenMP to accelerate memory transfer from 
% large array into small sub array and back
% 
% Matlab equivalent is : 
%   small_array = full_array(1:Npix_small(1),1:Npix_small(2),1:Npix_small(3));
%   full_array(1:Npix_small(1),1:Npix_small(2),1:Npix_small(3)) = full_array(1:Npix_small(1),1:Npix_small(2),1:Npix_small(3)) + small_array;


cd(fullfile( fileparts(mfilename('fullpath')), '..'))
addpath('tests')
addpath('utils')
addpath('./')
addpath(find_base_package)
utils.verbose(0)

%% if needed, recompile the mex functions manually
% cd cSAXS_matlab_base/+utils/private
% if verLessThan('matlab', '9.4')
%     mex -largeArrayDims 'CFLAGS="\$CFLAGS -fopenmp"' LDFLAGS="\$LDFLAGS -fopenmp" get_from_3D_projection.cpp
%     mex -largeArrayDims 'CFLAGS="\$CFLAGS -fopenmp"' LDFLAGS="\$LDFLAGS -fopenmp" add_to_3D_projection.cpp
% else
%     mex -R2018a 'CFLAGS="\$CFLAGS -fopenmp"' LDFLAGS="\$LDFLAGS -fopenmp" get_from_3D_projection_mex.cpp
%     mex -R2018a 'CFLAGS="\$CFLAGS -fopenmp"' LDFLAGS="\$LDFLAGS -fopenmp" add_to_3D_projection_mex.cpp
% end

Npix_full = [600,600,600]; 
Npix_small = [400,400,400]; 

utils.verbose(0,'--- get_from_3D_projection')
full_array = randn(Npix_full, 'single')+1i;
small_array = ones(Npix_small, 'like', single(1i));
positions = (10*rand(Npix_small(3),2));
indices = ([1:Npix_small(3)]);  % indices start from 1 

%% get 3D stack array "small_array" from 3D stack array "full_array" given the offsets "positions" and layers "indices"
utils.verbose(0,'Speed test MEX') 
for ii = 1:3
tic; utils.get_from_3D_projection(small_array,full_array,positions,indices); toc
end

assert(norm(small_array(:)-reshape(utils.get_from_3D_projection(small_array,full_array,positions,indices, false),[],1))==0, 'get_from_3D_projection MEX function does not provide exact results')


utils.verbose(0,'Speed test Matlab') 
for ii = 1:3
    tic; small_array = full_array(1:Npix_small(1),1:Npix_small(2),1:Npix_small(3)); toc
end

%% add 3D stack array "small_array" into 3D stack array "full_array" given the offsets "positions" and layers "indices"
utils.verbose(0,'--- add_to_3D_projection')
full_array = zeros(Npix_full, 'like', single(1i));
small_array = ones(Npix_small, 'like', single(1i));
positions = (10*rand(Npix_small(3),2));
indices = ([1:Npix_small(3)]);  % indices are starting from 1 !! 
add_values = true; 
utils.verbose(0,'=== add values atomic')
utils.verbose(0,'Speed test MEX') 
for ii = 1:3
    tic; utils.add_to_3D_projection(small_array,full_array,positions, indices,add_values, true);toc
end

assert(norm(full_array(:)-reshape(utils.add_to_3D_projection(3*small_array,zeros(Npix_full, 'like', single(1i)),positions, indices,add_values,true,false),[],1))==0, 'add_to_3D_projection MEX function does not provide exact results')

utils.verbose(0,'Speed test Matlab') 
tic; full_array(1:Npix_small(1),1:Npix_small(2),1:Npix_small(3)) = full_array(1:Npix_small(1),1:Npix_small(2),1:Npix_small(3)) + small_array; toc


%% set 3D stack array "small_array" into 3D stack array "full_array" given the offsets "positions" and layers "indices"
utils.verbose(0,'=== add values nonatomic')
add_values = true; 
utils.verbose(0,'Speed test MEX') 
full_array = zeros(Npix_full, 'like', single(1i));
for ii = 1:3
    tic; utils.add_to_3D_projection(small_array,full_array,positions, indices,add_values, false);toc
end
assert(norm(full_array(:)-reshape(utils.add_to_3D_projection(3*small_array,zeros(Npix_full, 'like', single(1i)),positions, indices,add_values,false,false),[],1))==0, 'add_to_3D_projection MEX function does not provide exact results')


utils.verbose(0,'Speed test Matlab') 
tic; full_array(1:Npix_small(1),1:Npix_small(2),1:Npix_small(3)) = full_array(1:Npix_small(1),1:Npix_small(2),1:Npix_small(3)) + small_array; toc
%% set 3D stack array "small_array" into 3D stack array "full_array" given the offsets "positions" and layers "indices"
utils.verbose(0,'=== set values')
add_values = false; 
utils.verbose(0,'Speed test MEX') 
full_array = zeros(Npix_full, 'like', single(1i));
for ii = 1:3
tic; utils.add_to_3D_projection(small_array,full_array,positions, indices,add_values);toc
end
assert(norm(full_array(:)-reshape(utils.add_to_3D_projection(small_array,zeros(Npix_full, 'like', single(1i)),positions, indices,add_values,false,false),[],1))==0, 'add_to_3D_projection MEX function does not provide exact results')

utils.verbose(0,'Speed test Matlab') 
tic; full_array(1:Npix_small(1),1:Npix_small(2),1:Npix_small(3)) = small_array; toc




%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2018 by Paul Scherrer Institute (http://www.psi.ch)    |
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

