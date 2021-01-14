%% UDPATE_MASK
% This small script guides you to update an alread existing mask for
% ptychography. The main tool for creating a new mask is
% beamline.create_mask, a GUI that lets you select bad/hot pixels. 
% UPDATE_MASK loads the data, specified by file_path, plots it and starts
% the GUI. Although you can create a 3D mask, i.e. a mask which varies from
% frame to frame, a 2D mask is sufficient for most datasets.

% You can load an already existing mask within the GUI.

close all


file_path = '~/Data10/eiger_4/S00000-00999/S00089/run_00089_000000000000.h5';
single_file = true;                     % if you have multiple files use * in file_path
H5Location = '/entry/data/eiger_4/';    % check the location within the h5 file in ptycho/+detector
orientation = [1 0 0];                  % check the detector orientation in ptycho/+detector



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the data

image_read_args = [];

image_read_args{1} = 'Orientation';
image_read_args{2} = orientation;
image_read_args{3} = 'OrientByExtension';
image_read_args{4} = false;

if ~single_file
    image_read_args{end+1} = 'IsFmask';
    image_read_args{end+1} = 1;
end

if ~isempty(H5Location)
    image_read_args{end+1} = 'H5Location';
    image_read_args{end+1} = H5Location;
end


data = io.image_read(file_path, image_read_args(:));

%% plot the data
figure(1),
plotting.imagesc3D(abs(log10(double(data.data)+1)));
colorbar
axis xy equal tight
colorbar 
title('Detector raw data')
colormap jet

%% iterative step for a mask update (add dead pixels to the current mask)
mask = beamline.create_mask;

%% check it again
figure (2),
imagesc(mask); axis equal tight xy
title('Final mask')

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

