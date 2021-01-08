% find_shift_3D_nonrigid - GPU accelerated weighted optical flow method 
%
% [shift,err] = find_shift_3D_nonrigid(vol_def, vol_ref, weight, downsample, smooth, regul)
%
% Inputs:
%    **vol_def     deformed volume        
%    **vol_ref     reference volume 
%    **weight      importance weights for each pixel
%    **downsample  downscale factor from the volume to DVF size 
%    **smooth      smoothness parameres for the recovered DVF
%    **regul       regularization preventing empty regions to have too large effect on the DVF estimate 
% Outputs: 
%    ++shift        calculated local shift for reference to match deformed volume  
%    ++err          error between reference and the deformed volume  

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


function [shift,err] = find_shift_3D_nonrigid(vol_def, vol_ref, weight, downsample, smooth, regul)
  
   import plotting.*
 
    
    % calculate error
    resid = vol_def-vol_ref;


    % apply high pass filtering
    resid = resid - utils.imgaussfilt3_fft(resid, 5);

    % calculate the error between the volumes 
    err = weight .* resid.^2;
    err = sqrt(mean(err(:))); 
    
    % avoid numerical instabilities
    weight = weight / mean(abs(resid(:)));

    Npix = size(vol_ref); 
    for i = 1:3
        ind_def{i} = gpuArray(linspace(1,Npix(i)/downsample, Npix(i))'); 
    end
    [X,Y,Z]= meshgrid(ind_def{:});

        
    for ax = 1:3
        % get gradient direction 
        vol_def_diff = math.get_img_grad_conv( vol_ref,2,ax);
        
        % estimate the optimal step 
        % GPU kernel merging 
        [num, denum]= arrayfun(@get_coefs,weight, resid, vol_def_diff);

        % bin the volume to make smoothing faster 
        num = utils.binning_3D(num, downsample); 
        denum = utils.binning_3D(denum, downsample);
        
        num =   padded_3D_smoothing(num, smooth/downsample/2); 
        denum = padded_3D_smoothing(denum, smooth/downsample/2); 

        % add some small regularization 
        denum = bsxfun(@plus, denum , regul*mean2(denum));
        
       
        shift{ax} = - num ./ denum;
        
        % run simple line search to refined the optimal step, ideal it should be close to 1
        shift_full = interp3(shift{ax}, X,Y,Z);
        update = shift_full.*vol_def_diff; 
       
        Nsteps = 10; 
        steps = logspace(0,1,Nsteps); 
        for ii = 1:Nsteps
            res = arrayfun(@get_residuum_err, weight, resid,update, steps(ii)); 
            err_tmp(ii) = gather(sum(sum(sum(res))));
            if ii > 1 && err_tmp(ii) > err_tmp(ii-1)
                break
            end
        end
        
        %% update the step         
        shift{ax}  = shift{ax} .* steps(math.argmin(err_tmp)); 
               

    end
    
    
end

function [num, denum]= get_coefs(W, resid, grad)
    % auxiliary function for fast GPU calculations 
    agrad = abs(grad); 
    W = W .* agrad; 
    % estimate the optimal step 
    num = W .* real(conj(resid) .* grad);
    denum = W .* agrad.^2;


end

function res = get_residuum_err(weight, resid, update, step)
    res = weight .* (resid + step.* update).^2;
end
            
function array = padded_3D_smoothing(array, smooth, split)
    % prevent periodic boundary issues for FFT conv smoothing 
    if nargin < 3 
        split = 1; 
    end
        
    Npad = ceil(min(size(array)/2, ceil(smooth/8)*16));
    array = padarray(array,[Npad(1),0,0],'symmetric','both');
    array = padarray(array,[0,Npad(2),0],'symmetric','both');
    array = padarray(array,[0,0,Npad(3)],'symmetric','both');
    
    array   = utils.imgaussfilt3_fft(array, smooth, split);

    array = array(Npad(1):end-Npad(1)-1, Npad(2):end-Npad(2)-1,Npad(3):end-Npad(3)-1); 

end

