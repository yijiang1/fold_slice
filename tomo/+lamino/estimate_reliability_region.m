% weight_sino = estimate_reliability_region(complex_projection, probe_size, subsample)
% Estimates the region where the complex projections are good enough to be unwrapped. Assumes that outside 
% of the measured FOV the amplitude of the object will be zero. This assumes that ptychography reconstruction
% did not add any values, due to those pixels never been reached by any probe element. Then the region is reduced
% by half the probe size using something akin to erosion. This function is only useful if the FOV is not square, 
% such as the case of the elliptical FOV of laminography.
% 
%  Pseudo code example: 
%  weights = imerode(abs(object) > 0, ones(probe_size/2))
% 
% Inputs:
%    **complex_projection       - (3D array) complex valued reconstructions 
%    **probe_size               - (int, int) size of the illumination probe 
%    **subsample                - (int) subsample the resulting array to save memory
% Outputs:
%   ++weight_sino                 - weights, 1 for full quality, 0<=W<1 for poor regions 


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


function weight_sino = estimate_reliability_region(complex_projection, probe_size, subsample)
    % simple reliability estimation based on amplitude of the
    % reconstruction 
    
    Npix = size(complex_projection); 
    Npix_new = ceil(Npix(1:2)/subsample/2)*2; 
    
    %% SOLVE THE PROBLEM IN LOW RESOLUTION 
    weight_sino = tomo.block_fun(@utils.interpolateFT, complex_projection,Npix_new, struct('use_GPU', false, 'use_fp16', false)); 
    probe_size = round(probe_size .* Npix_new ./ Npix(1:2)); 
    weight_sino =abs(weight_sino); 
    weight_sino = single(weight_sino > 0.1*quantile(weight_sino(:),0.9)); 
    
    %% only CPU is supported -> gather and move abck to GPU afterwards 
    weight_sino = gpuArray(utils.imcrop_outliers(gather(weight_sino)));     % leave only single largest compact object 

    
    [Y,X] = meshgrid(-ceil(probe_size(1)/2):floor(probe_size(1)/2), -ceil(probe_size(2)/2):floor(probe_size(2)/2));
    probe = (X/probe_size(1)*2).^2+(Y/probe_size(2)*2).^2 < 1; 
    
    kernel = probe/sum(probe(:)); 
    weight_sino = convn(weight_sino, kernel , 'same'); 
    weight_sino = weight_sino > 0.95; 
    weight_sino = uint8(imgaussfilt(single(weight_sino), 1)*255);    
    
end
