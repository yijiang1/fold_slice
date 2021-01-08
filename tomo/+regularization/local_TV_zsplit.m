% Vol = local_TV_zsplit(Vol, lambda, eps, Niter, use_chambolle)
% FUNCTION 
%    - split volume on smaller blocks for processing by local total variation 
%    - splitting is done automatically in order to fit to GPU memory 
% Inputs:
%     Vol - refined volume 
%     lambda - refinement step (tunning constant)
%     eps - regularization term (tunning constant)
%     Niter - number of iterations 
%     use_chambolle - if true use Chambolle method, if false use simple gradient descent solver 
% Recompile: 
%     mexcuda -output +regularization/private/local_TV_mex +regularization/private/TV_cuda_texture.cu +regularization/private/local_TV_mex.cpp

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



function Vol = local_TV_zsplit(Vol, lambda, eps, Niter, use_chambolle)
    
    keep_on_gpu = isa(Vol, 'gpuArray'); 

    gpu = gpuDevice;
    Nlayers = size(Vol,3); 
    Nblocks = ceil(numel(Vol)*4 / 1024e6);
    % empirically tested condition 
    Nblocks = max(Nblocks, ceil((numel(Vol)*4*14) / gpu.AvailableMemory));
    Nlayers_on_GPU = ceil(Nlayers/Nblocks); 
    Nblocks = ceil(Nlayers/Nlayers_on_GPU); 
    if Nblocks > 1
        for i = 1:Nblocks
            utils.progressbar(i, ceil(Nlayers/Nlayers_on_GPU))
            ind = (1+(i-1)*Nlayers_on_GPU):min(i*Nlayers_on_GPU, Nlayers);
            Vol_tmp = gpuArray(Vol(:,:,ind));  
            Vol_tmp = local_TV_mex(Vol_tmp,lambda,eps, Niter, use_chambolle); 
            if ~keep_on_gpu; Vol_tmp = gather(Vol_tmp);end
            Vol(:,:,ind) = Vol_tmp;
        end
    else
        Vol = local_TV_mex(Vol,lambda,eps, Niter, use_chambolle); 
        if ~keep_on_gpu; Vol = gather(Vol);end
    end
    
end
