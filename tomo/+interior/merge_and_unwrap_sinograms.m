% MERGE_WEIGHTS_INTERIORED_ARRAYS auxiliar function for interior tomography for block processing by tomo.block_fun 
% that adds up two arrays as X*W + (1-W)*Y, and unwrap resulting phase difference and
% return the phase 
%
%   phase_diff = merge_weights_interiored_arrays(stack_object, phase_diff_lowres, weights_interior,shift, ROI, param, exterior_weights_interior)
%
% Inputs: 
%   **stack_object        - complex valued projections 
%   **phase_diff_lowres   - phase gradient for the low resolution tomogram
%   **weights_interior              - relative weights_interiors of each region of the tomogram (denotes interior region ) 
%   **shift               -  shifts applied to the phase gradient 
%   **ROI                 - region to be unwrapped 
%   **param               - tomography parameter structure
%   **exterior_weights_interior     - scalar, relative weights_interior of the low res region during alignment
% Outputs: 
%   ++phase               - unwrapped phase of the low + high resolution sinogram together 
%   ++weights_interiors       - importance weights_interiors used for alignment, ie gives much smaller 
%                             weights_interior to the low resolution region compared to the interion projection part 
%                             based on the "exterior_weights_interior" input value 

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


function [phase_diff, weights_interiors] = merge_and_unwrap_sinograms(stack_object, tomo_lres, weights_interior,shift, theta, ROI, Nw_lres, param)

        

        Npix_lres = size(tomo_lres); 
        Npix_full = ceil(Nw_lres*param.resolution_ratio); 

        % get phase difference of the interior tomo 
        phase_diff_interior = math.get_phase_gradient_1D(stack_object,2,0); 
        
        avg_shift = round(mean(shift)); 
        Npix_partial = [size(stack_object,1), size(stack_object,2)]+2*ceil(max(abs(shift-avg_shift)));
        Npix_partial = min(Npix_partial, Npix_full);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% (1) UPSAMPLE WEIGHTS TO THE FULL RESOLUTION %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if isa(weights_interior, 'uint8') || ( isa(weights_interior, 'gpuArray') && strcmpi(classUnderlying(weights_interior),'uint8'))
            % convert to single precision 
            weights_interior = single(weights_interior-1)/255; 
        end
        if size(weights_interior,1) ~= size(stack_object,1) || size(weights_interior,2) ~= size(stack_object,2)
            % assume that the weights_interior are stored downsampled to same memory
            weights_interior = utils.interpolate_linear(weights_interior,  size(stack_object));
        end
       
        % pad arrays with enough space for save shifting 
        phase_diff_interior = utils.crop_pad(phase_diff_interior , Npix_partial);
        weights_interior    = utils.crop_pad(weights_interior , Npix_partial);

    
        
        % apply (shifts-avg_shift) on the interior tomogram and its weights_interiors
        [phase_diff_interior, weights_interior] = shift_interior_sinograms(phase_diff_interior, weights_interior, shift-avg_shift); 



        [cfg_lres, vectors_lres] = astra.ASTRA_initialize(Npix_lres,Nw_lres,theta, 90, 0, 1);

        % find optimal split of the dataset for given GPU 
        split = astra.ASTRA_find_optimal_split(cfg_lres);
        % forward projection model
        model = astra.Ax_partial(tomo_lres,cfg_lres, vectors_lres,split,'verbose', 0);

        % get phase difference
        phase_diff_lowres =  math.get_phase_gradient_1D( exp(1i*model), 2, 0);
   
        clear model
        % use FFT centred upsampling to keep it accurate    
        %phase_diff_lowres = utils.interpolateFT_centered(phase_diff_lowres / param.resolution_ratio, Npix_full,1); 
        
        % linear interpolation is much faster but it introduces bias
        bias = 1-1/param.resolution_ratio; 
        phase_diff_lowres = utils.imshift_fft(phase_diff_lowres,-[bias, bias]); 
        phase_diff_lowres = utils.interpolate_linear(phase_diff_lowres / param.resolution_ratio, Npix_full); 

        % pad arrays to full size and apply the average shift      
        phase_diff_interior = utils.imshift_fast(phase_diff_interior ,-avg_shift(1),-avg_shift(2), Npix_full); 
        weights_interior = utils.imshift_fast(weights_interior ,-avg_shift(1),-avg_shift(2), Npix_full);        

        alpha = 0.05; % avoid effects from weak regions of the weights 
        weights_interior = max(0, weights_interior-alpha)/(1-alpha); % remove artefacts from the FFT shift 

        % merge tomogram based on provided weights_interior 
        phase_diff = arrayfun(@merge_arrays,phase_diff_interior, phase_diff_lowres, weights_interior);
       
            
        if ~isempty(param.vert_range)
            % find optimal vertical range 
            Nlayers = length(ROI{1}); 
            vrange0 =  param.vert_range([1,end]); 
            vrange(2) = min(vrange0(2),Nlayers);
            vrange(1) = max(1,vrange0(1));

            % help with splitting in ASTRA (GPU memory limit)
            vrange_center = ceil(mean(vrange)); 
            estimated_split_factor = 4; 
            Nvert = floor((vrange(2)-vrange(1)+1) / estimated_split_factor)*estimated_split_factor;
            vrange(1) = ceil(vrange_center - Nvert/2); 
            vrange(2) = floor(vrange_center + Nvert/2-1); 
            vrange = vrange(1):vrange(2);
        else
            vrange = ':';
        end

    

        % select only object_ROI -> make it easily splitable to GPU and avoid
        % edge artefacts  
        phase_diff       = phase_diff(ROI{1}(vrange),ROI{2},:);
        weights_interior = weights_interior(ROI{1}(vrange),ROI{2},:);


        % apply binning 
        if param.binning > 1
            phase_diff = utils.interpolateFT_centered(phase_diff,ceil(size(phase_diff)/param.binning/2)*2, 1); % accurate interpolation using FFT
        end
        
        if strcmpi(param.unwrap_data_method, 'fft_1d')
            phase_diff = math.unwrap2D_fft(phase_diff,2,param.air_gap,0);
        end

        % provide weights_interiors to be used for alignment 
        weights_interiors = param.exterior_weight + (1-param.exterior_weight)*weights_interior; 
        
        
        weights_interiors = utils.interpolate_linear(weights_interiors,  ceil(size(weights_interiors)/10));

        % keep in uint8 to save memory 
        weights_interiors = uint8(weights_interiors*255); 
        
        
        
end

function C = merge_arrays(A, B, W)
    W = max(0, min(1,W)); 
    C = A.*W + (1-W).*B;
end


% FUNCTION [phase_diff_interior,weights_interior ] = shift_interior_sinograms(phase_diff_interior, weights_interior, shift )
% auxiliar function for interior tomography 
% Inputs: 
%   phase_diff_interior - phase gradient for the interior tomo 
%   weights_interior  - relative weights_interiors of each region of the tomogram (denotes interior region ) 
%   shift -  shifts applied to the phase gradient 
% Outputs: 
%   phase_diff_interior - phase difference for only the interior part
%   weights_interior - reliability weights_interiors estimated from the original field of view and shifted accordingly

function [phase_diff_interior,weights_interior ] = shift_interior_sinograms(phase_diff_interior, weights_interior, shift )
    
    % shift phase and weights_interiors in parallel to save time and avoid numerical problems at the edges 
    phase_diff_interior = utils.imshift_fft(weights_interior .* exp(1i*phase_diff_interior),shift);
        
    % get amplitude (weights_interiors)
    weights_interior = abs(phase_diff_interior);

    % get phase = phase diff 
    phase_diff_interior = angle(phase_diff_interior);

end
