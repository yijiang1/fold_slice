% ALIGN_TOMO_XCORR Cross-correlation alignment 
% [total_shift, variation, variation_aligned]  = align_tomo_Xcorr(object, angles, par, varargin)
%
% Inputs:
%   **object    - complex valued projections to be aligned
%   **angles    - corresponding angles (used only for sorting the projections)
%   **par       - parameter structure with: 
% *optional*   (use value in param as default)
%   **filter_pos - parameter for high pass filtering of the shifts
%   **filter_data - highpass filtering of the data 
%   **max_iter  - maximal number of iterations 
%   **binning   - binning applied on the data 
%   **ROI       - region of interest 
% *returns*
%   ++total_shift   - optimal shift of the projections 
%   ++variation     - local variation of the measured data 
%   ++variation_aligned  - aligned variation of the measured data 


%        
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

function [total_shift, variation, variation_aligned]  = align_tomo_Xcorr(object_0, angles, par, varargin)

    import utils.*
    import math.*
    utils.verbose(struct('prefix', 'align'))
    
    parser = inputParser;
    parser.addParameter('binning',  1 , @isnumeric )
    parser.addParameter('max_iter', 1 , @isnumeric )
    parser.addParameter('filter_pos', 50 , @isnumeric )
    parser.addParameter('filter_data', 0.05 , @isnumeric )
    parser.addParameter('ROI', {} , @iscell )


    parser.parse(varargin{:})
    r = parser.Results;

    % load all to the param structure 
    for name = fieldnames(r)'
        if ~isfield(par, name{1}) || ~ismember(name, parser.UsingDefaults) % prefer values in param structure 
            par.(name{1}) = r.(name{1});
        end
    end
    
    if isreal(object_0)
        error('Complex object expected')
    end
    % binning to speed up the calculation, anyway we need only low
    % resolution guess 
    
    weights = par.illum_sum ./ (par.illum_sum+1e-1*max(par.illum_sum(:))); 
    variation = tomo.block_fun(@get_variation_field, object_0,par.binning, weights(par.ROI{:}), struct('use_GPU', true, 'ROI', {par.ROI}, 'use_fp16', false));     
    variation = real(variation);  % the real() function must be applied outside of get_variation_field to get betetr performance, it seems like a bug in matlab 
    
    [~,ind_sort] = sort(angles);
    [~,ind_sort_inv] = sort(ind_sort);

    
    % get some initial guess of the relative shifts 
    [Nx, Ny, Nangles] = size(variation);
    total_shift = zeros(length(angles),2);

    
    if ~par.is_laminography
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find center of rotation from comparison between 0 and 180 deg frame %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ind = math.argmin(abs(angles(1)+180 - angles)); % find the opposite angle to 0 degrees
        avg_step = 5*median(diff(sort(angles))); 
        if abs(angles(ind)-(angles(1)+180)) < avg_step
            obj_tmp = variation(:,:,[1,ind]); 
            obj_tmp(:,:,2) = fliplr(obj_tmp(:,:,2));
            fvar = filtered_FFT(obj_tmp,[0,0],par); 
            CoR_offset = 0.5 * find_shift_fast_2D(fvar(:,:,1),fvar(:,:,2),0, false); 
            % apply the estimated CoR offset to the total shift 
            total_shift(:,1) = total_shift(:,1) - CoR_offset(1); 
        else
           utils.verbose('Missing mirror projection for initial projection, skipping CoR estimation') 
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % apply crosscorrelation between subsequent frames to remove the relative shifts %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for iter = 1:par.max_iter
        utils.progressbar(iter, par.max_iter);  
        

        % apply the estimated shift to the data to get better next step  
        % AND
        % use local variation in order to suppress issues
        % caused by low spatial freq. errors in the ptychography
        % reconstruction
        fvar = tomo.block_fun(@filtered_FFT,variation,total_shift, par, struct('verbose_level', 0));     
  
        % compare subsequent slices 
        frame_ref = fvar(:,:,ind_sort);
        frame_align = fvar(:,:,circshift(ind_sort,-1)); 

        % align the first frame with flipped version of the last frame 
        if  ~par.is_laminography && abs(mod(angles(1) - angles(Nangles)-180, 360)) < 2*median(abs(diff(angles)))
            frame_ref(:,:,end) = fft2(fliplr(ifft2(frame_ref(:,:,end))));
        end
        
        % find the optimal relative shift between projections and the adjanced angles 
        relative_shifts = tomo.block_fun(@find_shift_fast_2D,frame_ref,frame_align,0, false,struct('verbose_level', 0));     
        relative_shifts = circshift(relative_shifts,1);        
        

        % avoid too fast jumps by limiting the maximal step size per iteration 
        max_shift = max(10 ,3*mad(relative_shifts)); 
        relative_shifts = min(max_shift, abs(relative_shifts)) .* sign(relative_shifts);

        % long term drifts cannot be trusted 
        cum_shift = cumsum(relative_shifts);
        cum_shift = cum_shift - mean(cum_shift);

        % minimize the shift amplitude that is needed 
        % long term drifts cannot be trusted 
        if ~isinf(par.filter_pos)
            for i = 1:2
                smooth = ceil(par.filter_pos/2)*2+1; 
                % get properly smoothed cumulative shift 
                smoothed_shift = conv(cum_shift(:,i),ones(smooth,1)/smooth, 'same') ./ conv(ones(Nangles,1),ones(smooth,1)/smooth, 'same'); 
                % subtract it from cum_shift
                cum_shift(:,i) = cum_shift(:,i) - smoothed_shift; 
            end
        end


        total_shift = total_shift + cum_shift(ind_sort_inv,:); 
        % limit the maximal shift to 3* mean absolute devition of all the
        % position -> prevent outliers 
        total_shift = min(6*mad(total_shift,1), abs(total_shift)) .* sign(total_shift); 
        
      
        % draw position shifts 
        plotting.smart_figure(12)
        clf()
        plot(par.scanstomo, (total_shift(ind_sort,:)-mean(total_shift)) * par.binning,'.')
        legend({'Horizontal', 'Vertical'})
        xlabel('Scan number')
        axis tight 
        grid on
        title('Total shift estimated by initial cross-correlation')

        drawnow
        
        
        if 3*mad(cum_shift(:)) < par.precision && ~debug()
            break
        end
    end
    
    % return shifted projections for user to judge the alignment quality 
    variation_aligned =  utils.imshift_fft(variation, total_shift); 

    
    total_shift = total_shift .* par.binning; 
    
    if ~debug()
        total_shift = round(total_shift); % the accuracy is too low for subpixel shift
    end
    
end

function variation = get_variation_field(object, binning, weights)
    % just auxiliar function to move computations on GPU 

    import utils.*
    import math.*
    
    dX = convn(object,[-1,1], 'same'); 
    dY = convn(object,[-1,1]', 'same'); 

    if isa(dX, 'gpuArray')
        variation = arrayfun(@aux_fun,dX, dY, object);
    else
        variation = aux_fun(dX, dY, object);
    end

    % remove phase ramp artefacts 
    variation([1,end],:,:) = variation([2,end-1],:,:);    
    variation(:,[1,end],:) = variation(:,[2,end-1],:) ; 
    
    % crop values exceeding limits, important mainly for laminography where
    % the field of view can contain weakly illuminated (ie very noisy) regions
    mean_variation = mean2(variation .* weights) ./ mean2(weights); 
    dev_variation = sqrt(mean2((variation-mean_variation).^2 .* weights) ./ mean2(weights)); 
    
    variation = min(variation, mean_variation + 1*dev_variation); 
      
    
    % decimate data to lower resolution, smoothing needs to be applied before downsampling
    % smoothed array works better than binning
    variation = utils.imgaussfilt3_conv(variation, [2*binning,2*binning,0]);
    boundary_correction = utils.imgaussfilt3_conv(ones(size(object,1), size(object,2), 'like', object), [2*binning,2*binning,0]); 
    variation = variation ./ boundary_correction; 
    variation = variation(1:binning:end,1:binning:end,:); 
    
      
end

function variation = aux_fun(dX, dY, object)

    % get total variation -> better results in alignment than raw
    % object
    variation = sqrt(abs(dX).^2 + abs(dY).^2); 
    
    % ignore regions with low amplitude 
    variation = variation .* abs(object); 
        
end


function img = filtered_FFT(img, shift, par)

    [nx, ny, ~] = size(img);

    img = utils.imshift_fft(img ,shift);
            
    % suppress edge effects of the registration procedure
    spatial_filter = tukeywin(nx,0.3) * tukeywin(ny,0.3)'; 
    img = img - mean(img(:)); 
    img = img .* spatial_filter;

    % precalculate FFT
    img = fft2(img);

    % remove low frequencies (e.g. phase ramp issues) 
    if par.filter_data > 0 
        [X,Y] = meshgrid( (-nx/2:nx/2-1), (-ny/2:ny/2-1));
        spectral_filter = fftshift(exp(1./(-(X.^2+Y.^2)/(mean([nx,ny])*par.filter_data)^2)))';
        img =  img.* spectral_filter;
    end
end

