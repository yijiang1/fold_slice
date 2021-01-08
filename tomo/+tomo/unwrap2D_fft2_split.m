%  UNWRAP2D_FFT2_SPLIT simple and very fast 2D phase unwrapping with autosplitting for GPU
%  It applies iterativelly utils.unwrap2D_fft2 and enforces constrains by
%  remove_sinogram_ramp, if abs(angle(img .* exp(-1i*phase))) < 2  , use
%  phase = phase + angle(img .* exp(-1i*phase)) for exact unwrapping 
%
%  Method: estimate phase gradients dX, dY, and perform 2D complex
%  integration as for DIC method to get phase (as in p = phase_from_dpc(dpcx,dpcy,varargin) function)
%  
%  method is similar (but not identical) to 
%  Sam Jeught, Jan Sijbers, and Joris Dirckx. "Fast Fourier-based phase unwrapping on the graphics processing unit in real-time imaging applications." Journal of Imaging 1.1 (2015): 31-44.
%
%  [varargout] = unwrap2D_fft2_split(img,  empty_region, polyfit_order, weights,  GPU_list, ROI, Niter)
%
% Inputs: 
%   **img                -  either complex valued image or real valued phase gradient 
%   **empty_region       - 2x1 or 1x1 vector, size of empty region assumed around edges for phase offset removal , default = []
%   **polyfit_order      -   -1 = dont assume anything about the removed phase,
%                       subtract linear fit a*x+b for each horizontal line in order to satisfy
%                       that values in the empty_region are zero 
%                      0 = (default) assume that it is constant offset and minimize values in the empty_region
%                      1 = assume phase ramp it is 2D plane. monimize values in empty_region
%   **weights            - reliability weights from 0 to 1 ( default = 1), can be just a function
%                         handle taking as input "img" array or a downsampled array that will 
%                        be fourier interpolated before unwrapping, 
%   **GPU_list           - list of used GPUs, default = current GPU 
%   **ROI                - unwrapped region, default ROI = {':',':'};  
%   **preprocess_fun     - apply custom function on "img" before processing 
%   **Niter              - maximal number of unwrapping refinement interations, default = 5
% *returns*
%   ++phase              - unwrapped phase
%  
% Examples: 
%    x = linspace(0, 10, 100);
%    xc = exp(10*sin(x).*cos(x')); % make some 2D complex valued array 
%    xc = repmat(xc, 1,1, 100) ; % just show that it works for stacked inputs 
%    x_unwrapped = unwrap2D_fft2_split(img);   % simplest case, no boundary conditions are applied


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




function [varargout] = unwrap2D_fft2_split(img,  empty_region, polyfit_order, weights_0, GPU_list, ROI, preprocess_fun , Niter)

    import utils.* 
    import math.* 

    if isreal(img) && ~isa(img, 'uint32')
        error('Complex-valued input array was expected')
    end

    if nargin < 3
        polyfit_order = 1; 
    end
    if nargin < 2
        empty_region = [];
    end
    if nargin < 4
        weights_0 = 1;
    end
    if nargin < 5
        GPU_list = []; % use default GPU 
    end
    if nargin < 6 || isempty(ROI)
        ROI = {':',':'};  % unwrap only a small ROI from the full complex array
    end

    if nargin < 7
        preprocess_fun = [];
    end
    if nargin < 8
        Niter = [] ;
    end
    gpu = gpuDevice; 
    if ~isempty(GPU_list) && ~ismember(gpu.Index, GPU_list)
        gpu = gpuDevice(GPU_list(1)); 
    end


    [Nx,Ny] = size(img(ROI{:},1)); 
    Nz = size(img,3); 
    
    if gpuDeviceCount
        gpu = gpuDevice; 
        if ~ismember(gpu.Index, GPU_list) && ~isempty(GPU_list)
            if isa(img, 'gpuArray')
                error('Non gpuArray input expected, change of GPU id will reset GPU memory content')
            end

            gpu = gpuDevice(GPU_list(1)); 
        end
        AvailableMemory = gpu.AvailableMemory; 
    else
        % run in RAM 
        AvailableMemory = utils.check_available_memory*1e6; 
    end
    
    Nblocks =  ceil( (2e9+  10 *8* (Nx+128)*(Ny+128)*size(img,3))  /  AvailableMemory) ; 
    Nblocks = max(Nblocks,  (Nx+64)*(Ny+64)*size(img,3)  / double(intmax('int32')));
    % avoid issues with rouding of Nz
    Nblocks = ceil(Nz / floor(Nz/Nblocks)); 
    

    if ~isempty(weights_0) && isnumeric(weights_0) && (ismatrix(weights_0)  || any(size(img) ~= size(weights_0)))
       if any([ size(img,1),size(img,2)] ~= [size(weights_0,1),size(weights_0,2)])
          for i = 1:2
                wROI{i} = unique(ceil(ROI{i}*size(weights_0,i) / size(img,i)));
          end
       else
           wROI = ROI; 
       end
       weights_0 = weights_0(wROI{:},:); 
    end
    
    params = struct('Nblocks', Nblocks, 'GPU_list', GPU_list, 'ROI', {ROI}, 'use_GPU', gpuDeviceCount > 0, 'use_fp16', false, 'move_to_GPU', false); 
    
    varargout = cell(nargout,1); 
    [varargout{:}] = tomo.block_fun(@unwrap2D_fft2_worker, img, empty_region,weights_0,polyfit_order,preprocess_fun, Niter, params);
        
end


function [phase_block, residues_block] = unwrap2D_fft2_worker(img_block, empty_region,weights_0,polyfit_order,preprocess_fun, Niter)
    import utils.*
    import math.*
    Npix = size(img_block); 
    if isempty(weights_0)  || isscalar(weights_0)
        weights = ones(size(img_block,1), size(img_block,2), 'single'); 
    elseif  isa(weights_0, 'function_handle')
        weights = weights_0(img_block);
    elseif isnumeric(weights_0) && any(Npix(1:2) ~= [size(weights_0,1),size(weights_0,2)])
        weights_0 = gpuArray(weights_0); 
        weights_0 = single(weights_0) / 255; 
        weights = utils.interpolate_linear(weights_0, Npix(1:2)); 
    else
        weights = weights_0; 
    end
    weights = Garray(weights); 
    img_block = Garray(img_block); 
    if any(~isfinite(img_block(:)))
        error('Unwrapped complex array contains nan/inf values')
    end
    % apply custom function if provided  
    if ~isempty(preprocess_fun)
        img_block = preprocess_fun(img_block); 
    end
    weights = max(0,weights) / max(weights(:)); 

    
    %phase_block = unwrap2D_fft2(img_block, empty_region,0,weights,polyfit_order);
    
    % find residua, it is computationally cheap 
    residues_block = abs(findresidues(img_block)) .* weights(2:end,2:end,:) > 0.1; 
    residues_block = uint8(residues_block);  % add_to_projection MEX function does not support logicals -> use uint8 which has the same size in matlab 
    
    % decide how many refinement iterations 
    if isempty(Niter)
        if any(residues_block(:))
           %% internal variable to set number of iterative refinements 
           Niter = 10; 
        else
           Niter = 5; 
        end
    end

    % initialize resulting phase 
    phase_block = 0; 

    % perform several iterations to refine the quality
    W = weights; 
    for iter = 1:Niter
        if iter == 1 
            % initial unwrapping 
            img_block_resid =img_block; 
        else
            img_block_resid =img_block.*exp(-1i*phase_block); 
        end
        %% FOR DEBUGGING 
%         plotting.imagesc3D( W.*angle(img_block.*exp(-1i*phase_block)), 'init_frame', 1)
%         axis xy
%         colormap hsv(1024)
%         colorbar 
%         caxis([-pi,pi])
%         title(sprintf('Iter %i', iter))
%         pause(1)
%         drawnow 

        % check that unwrapping is really needed 
        [a_resid,~,~,~,c_factor] = utils.stabilize_phase(img_block_resid, 'fourier_guess', false, 'weight', W); 
        a_resid = angle(a_resid); 
        
        if all(all(all(abs(W.* (a_resid) )< 2 )))
            % if the data are nice, make !! exact !! unwrapping and finish 
            phase_block = phase_block + W.*(a_resid-c_factor);
            if ~isempty(empty_region)
                % but still be sure to properly remove phase ramp / offset 
                phase_block = remove_sinogram_ramp(phase_block,empty_region, polyfit_order);  
            end
            return
        end
        clear a_resid
        
        phase_block = phase_block + unwrap2D_fft2(img_block_resid,[],0, W, polyfit_order);
        
        if ~isempty(empty_region)
            % but still be sure to properly remove phase ramp / offset 
            phase_block = remove_sinogram_ramp(phase_block,empty_region, polyfit_order);  
        end

    end     
end

