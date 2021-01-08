% ALIGN_TOMO_CONSISTENCY_XCORR cross-correlation of measured projections with the projected model 
% is it much more robust to outliers than the optimization based methods, 
% ! but in some cases it can make the reconstruction worse ! 
% also only model projections are shifted in contrast to align_tomo_consistency_linear where the measured projectins are shifted
% issues that majority of the points are correct
%
% xcorr_shift_total = align_tomo_consistency_Xcorr(sinogram_0, sino_weights, angles, shift_0, binning, Npix, par, varargin)
%
% Inputs:
%   **sinogram_0        - unwrapped phase or phase difference 
%   **angles            - corresponding angles (used only for sorting the projections)
%   **shift_0           - initial guess of the sinogram shifts 
%   **binning           - binning used for calculations to speed it up 
%   **Npix              - size of the tomogram   
%   **par               - tomography parameter structure , 
%   **varargin          - list of all paremters is described in code 
% *returns*: 
%   ++xcorr_shift_total - additional shift that nees to be added to shift_0 input in order to maximize this self-consitency alignment method 


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



function xcorr_shift_total = align_tomo_consistency_Xcorr(sinogram_0, sino_weights, angles, shift_0, binning, Npix, par, varargin)


    import tomo.*
    import utils.*
    import math.*
    import plotting.*
    utils.verbose(struct('prefix', 'align'))

    
    parser = inputParser;
    parser.addParameter('align_vertical',  true , @islogical )
    parser.addParameter('align_horizontal', false , @isnumeric )
    parser.addParameter('apply_positivity',  true , @islogical )  % apply_positivity , useful for tomography, should be avoided for laminography 
    parser.addParameter('high_pass_filter',  0.02 , @isnumeric )  % high pass filter applied on the xcorrelated model and data to remove effect of low spatial freq. errors, ie phase ramp , smaller value == less filtering 
    parser.addParameter('binning',  4 , @isint )                   % binning applied on the projections to speed up shift estimation
    parser.addParameter('lamino_angle',  90 , @isnumeric )         % laminography title (with respect to the beam ))
    parser.addParameter('tilt_angle',  0 , @isnumeric )            % rotation of camera around axis of the beam 
    parser.addParameter('skewness_angle',  0 , @isnumeric )        % skewness of camera around axis of the beam 
    parser.addParameter('pixel_scale',  [1,1] , @isnumeric )       % pixel scale along each axis 
    parser.addParameter('unwrap_data_method',  'fft_1D' , @isstr ) % assume  that inputs is phase derivative, accepted values: 'fft_1d' for input phase difference, 'none' for already unwrapped input
    parser.addParameter('verbose',  1 , @isnumeric )               % change verbosity of the code 
    parser.addParameter('filter_type',  'ram-lak' , @isstr )       % FBP settings 
    parser.addParameter('freq_scale',  '1' , @isnumeric )          % FBP settings  
    parser.addParameter('plot_results',  true , @islogical )       % plot results 
    parser.addParameter('is_laminography',  false , @isnumeric )   % change default behavoir for laminography reconstructions
    parser.addParameter('CoR_offset',  [] , @isnumeric )           % offset of the center of rotation , empty == auto  
    parser.addParameter('Niter',  5 , @isnumeric )                 % number of ierations of the Xcorr refinement 

    parser.addParameter('selected_roi',  {} , @iscell ) % field of view considered for alignment for laminography
    parser.addParameter('vert_range',  [] , @isnumeric) % vertical range considered for alignment for standard tomo 

    
    parser.parse(varargin{:})
    r = parser.Results;

    % load all to the param structure 
    for name = fieldnames(r)'
        if ~isfield(par, name{1}) || ~ismember(name, parser.UsingDefaults) % prefer values in param structure 
            par.(name{1}) = r.(name{1});
        end
    end
    
    assert(isreal(sinogram_0), 'Input projections has to be real-valued, either unwrapped phase or phase difference')

    
    verbose(1,'==== Starting self-consistent cross correlation  ==== ')
        
    relax = 0.8; 
    
    if ~isempty(par.vert_range) && isempty(par.selected_roi)
        vrange = ceil(par.vert_range([1,end])/16)*16 + [0,-1];  % help with splitting in ASTRA (GPU memory)
        if isempty(vrange(1):vrange(2)) 
            error('Too small par.vert_range, extend the alignment range'); 
        end 
        if length(vrange(1):vrange(2)) < 10*binning && par.align_vertical
            error('Too small par.vert_range for vertical alignment, extend the alignment range'); 
        end 
        par.selected_roi = {vrange(1):vrange(2),':'};
    end
    [Nlayers,width_sinogram,Nangles] = size(sinogram_0);

    verbose(1,['Binning: ', num2str(binning)])
    sinogram_0 = tomo.block_fun(@imreduce,sinogram_0, par.selected_roi, binning, [Nlayers,width_sinogram]);
    if ~isempty(sino_weights) 
        if isa(sino_weights, 'uint8')
            sino_weights = single(sino_weights)/255; % load the weights from uint8 format
        end
        assert(all(sum(sum(sino_weights)) > 0), sprintf('Provided "weights" contain %i projections with empty mask', sum(sum(sum(sino_weights))==0)))
        sino_weights = tomo.block_fun(@imreduce,sino_weights, par.selected_roi, binning, [Nlayers,width_sinogram], struct('full_block_size', [Nlayers,width_sinogram,Nangles]));
    else
        sino_weights = 1 ;
    end
    [Nlayers,width_sinogram,Nangles] = size(sinogram_0);
    
    if gpuDeviceCount
        gpu = gpuDevice;
        if isempty(par.GPU_list)
            par.GPU_list = gpu.Index;
        end
        GPU_list =  par.GPU_list(1); 
        gpu_split = 1; 

        % split it among multiple GPU only for larger datasets (>1000MB)
        if numel(sinogram_0) * 4 > 1000e6 && length(par.GPU_list) > 1      
            gpu_split =  max(1,length(par.GPU_list)); 
            GPU_list = par.GPU_list; 
        end
        if gpu.Index ~= GPU_list(1)
            gpuDevice(GPU_list(1));
        end
    else
        gpu_split = 1; 
        GPU_list = []; 
    end
    
    
    
    Npix = ceil(Npix/binning);
    if isscalar(Npix)
        Npix = [Npix, Npix, Nlayers];
    elseif length(Npix) == 2
        Npix = [Npix, Nlayers];
    end
    shift_0 = shift_0 / binning;
    
    
    %% unwrap complex data 
    if ~strcmpi( par.unwrap_data_method, 'none')
        verbose(1,'Unwrapping')

        sinogram_unwrapped = tomo.block_fun(@unwrap_data, sinogram_0,  par.unwrap_data_method, par.air_gap/par.binning); 
    else
        sinogram_unwrapped = sinogram_0; 
    end
            
    %% get configuration for ASTRA
    % !! important for binning => account for additional shift of the center
    % of rotation after binning, for binning == 1 the correction is zero
    if strcmpi(par.unwrap_data_method, 'none')
         % direclty unwrapped phase 
         rotation_center = [Nlayers, width_sinogram]/2 + (1-1/par.binning);
    else  % phase difference 
         rotation_center = [Nlayers, width_sinogram]/2 - (1-1/par.binning);
    end
    if ~isempty(par.CoR_offset) 
        rotation_center = rotation_center  + par.CoR_offset/par.binning; 
    end

    [cfg, vectors_0] = ...
        astra.ASTRA_initialize(Npix,[Nlayers, width_sinogram],angles,par.lamino_angle,par.tilt_angle,1, rotation_center); 
    %% find optimal split of the dataset for given GPU 
    split = astra.ASTRA_find_optimal_split(cfg, gpu_split);
  
    vectors = vectors_0; 
    xcorr_shift_total = zeros(Nangles,2); 

    %% choose solver
    if gpuDeviceCount == 0 
        tomo_solver = @FBP_CPU; 
        padding = []; 
    elseif ~par.is_laminography
        % use simple z-axis splitting 
        tomo_solver = @FBP_zsplit ; 
        padding = 0; 
    else
        % use more general solver, can be less efficient 
        tomo_solver = @FBP; 
        padding = 'symmetric'; 
    end

    win = tukeywin(width_sinogram, 0.2)';  % avoid edge issues 
    if Nlayers > 10 && (par.align_vertical ); win = tukeywin(Nlayers, 0.2)*win; end 
    sino_weights = sino_weights .* win; 

    
    for ii = 1:par.Niter 
        verbose(1,'Iter %i/%i ', ii, par.Niter)

    
        %% apply initial shift to the geometry 
        vectors(:,4:6) = vectors_0(:,4:6) + ...
            bsxfun(@times,vectors_0(:,10:12), shift_0(:,2) + xcorr_shift_total(:,2)) +...
            bsxfun(@times,vectors_0(:,7:9),   shift_0(:,1) + xcorr_shift_total(:,1));

        %% FBP method 
        rec  = tomo_solver(sinogram_unwrapped, cfg, vectors, [1,1,gpu_split],...
            'valid_angles',par.valid_angles, ...
            'GPU', GPU_list, 'split_sub', split, 'verbose', 0,...
                'filter', par.filter_type, 'filter_value', par.freq_scale, 'padding', padding);


        if par.apply_positivity
            rec = max(0, rec); 
        end
                

        % forward projection 
        if gpuDeviceCount
            sinogram_corr = Ax_sup_partial(rec,cfg, vectors,[1,1,max(1,length(par.GPU_list))],...
                'GPU', GPU_list, 'split_sub', split, 'verbose', 0);
        else
            sinogram_corr = radon_wrapper(rec, cfg, vectors);
        end
        % filter data to remove low freq. errors 
        
        if strcmpi(par.unwrap_data_method, 'fft_1d')
            % fft_1d was used for unwrap -> input sinogram is phase_difference
            sinogram_corr = - tomo.block_fun(@math.get_phase_gradient_1D,sinogram_corr);
            % reduce effect of too sharp features 
            sinogram_corr = min(0.1, abs(sinogram_corr)) .* sign(sinogram_corr); 
            sinogram = min(0.1, abs(sinogram_0)) .* sign(sinogram_0); 
        elseif strcmpi(par.unwrap_data_method, 'none')
            % suppress lowest spatial frequencies such as phase ramp from effecting alignment 
            sinogram_corr = sinogram_corr - tomo.block_fun(@utils.imgaussfilt2_fft, sinogram_corr, width_sinogram / 10); 
            sinogram = sinogram_unwrapped - tomo.block_fun(@utils.imgaussfilt2_fft, sinogram_unwrapped, width_sinogram / 10); 
        else
            error('Unsupported unwrap_data_method')
        end
        

        % perform 2D cross-correlation to find the shift 
        xcorr_shift = utils.find_shift_fast_2D(sinogram_corr.*sino_weights,sinogram.*sino_weights, par.high_pass_filter, 'full_range');
        xcorr_shift = gather(xcorr_shift); 
        % remove constant offsets 
        xcorr_shift = xcorr_shift - median(xcorr_shift);
        
        % avoid too large jumps 
        xcorr_shift = xcorr_shift * relax;        
        xcorr_shift = min([width_sinogram,Nlayers] / 4, abs(xcorr_shift)) .* sign(xcorr_shift); 
        
        if any(isnan(xcorr_shift(:)))
            warning('Crosscorrelation alignment resulted in NaN positions')
            keyboard
        end

        if ~par.align_vertical
            xcorr_shift(:,2) = 0;
        end
        if ~par.align_horizontal
            xcorr_shift(:,1) = 0;
        end
        
        
        xcorr_shift_total = xcorr_shift_total + gather(xcorr_shift); 
        
        if any(xcorr_shift_total(:)~=0) && par.plot_results
            [angles_sort, ang_sort] = sort(angles);
            plotting.smart_figure(6464)
            clf()
            warning('off', 'MATLAB:legend:IgnoringExtraEntries')
            if par.align_horizontal
                subplot(1,2,1)
                plot(angles_sort, xcorr_shift_total(ang_sort,1)*binning, '.r')
                subplot(1,2,2)
                plot(angles_sort, shift_0(ang_sort,1)*binning, '.r')
            end
            if par.align_vertical
                subplot(1,2,1)
                hold on 
                plot(angles_sort, xcorr_shift_total(ang_sort,2)*binning, '.b')
                hold off
                subplot(1,2,2)
                hold on 
                plot(angles_sort, shift_0(ang_sort,2)*binning, '.b')
                hold off
            end
            subplot(1,2,1)
            legend({'Horizontal ', 'Vertical'});
            title('Xcorr update shift')
            xlim([min(angles), max(angles)])
            xlabel('Angles')
            ylabel('Shift [px]')
            subplot(1,2,2)
            legend({'Horizontal ', 'Vertical'});
            title('Initial shifts')
            xlim([min(angles), max(angles)])
            xlabel('Angles')
            ylabel('Shift [px]')
            
            suptitle('Consitency-based cross-correlation position correction')
            
            warning('on', 'MATLAB:legend:IgnoringExtraEntries')

            drawnow 
        end
        
        if max(abs(xcorr_shift))*binning < 1
            % convergence reached 
            progressbar(par.Niter, par.Niter)
            break 
        end


    end
    
    xcorr_shift_total = xcorr_shift_total * binning;


end

function img = imreduce(img, ROI, binning, Npix)
    % auxiliary function that first upsamples the provided array to Npix
    % size, then crops it down to ROI region and finally perform
    % downsampling by interpolateFT that provides more precise results than
    % binning
    % !! merging all these operations into one function makes it much more efficient when processed on GPU via tomo.block_fun 
    %
    % Inputs
    % ** img   - input image block 
    % ** ROI   - cells of indices for selected ROI , ie. {10:500, 63:300}, ROI cropping is applied after first upsampling to Npix
    % ** binning - integer or 2x1 vector, binning the img 
    % ** Npix - target size of img before binning - useful to save mask in low resolution and interpolate it when needed 
    % returns: 
    % ++img  - processed image block 
    
    import math.*
    import utils.*
    
    isReal = isreal(img); 
    
    % if needed upsample to the size of the projection 
    img = utils.interpolate_linear(img, Npix); 

    % crop the FOV after shift and before "binning"
    if ~isempty(ROI)
       img = img(ROI{:},:);   % crop to smaller ROI if provided  
                              % apply crop aftetr imshift_fft
    end
    
    Np = size(img); 
    % perform FT interpolation instead of binning 
    img = interpolateFT(img, ceil(Np(1:2)/binning/2)*2);
    if isReal; img = real(img); end
end


function sinogram = unwrap_data(sinogram, method, boundary)
    % auxiliary function to perform data unwrapping, see
    % math.unwrap2D_fft for detailed help 
    %  ** sinogram  - unwrapped arrays 
    %  ** method  - none or fft_1d unwrap method 
    %  ** boundary - air region for zero boundary condtion 
    
    switch lower(method)
        case 'fft_1d'
            % unwrap the data by fft along slices 
            sinogram = -math.unwrap2D_fft(sinogram, 2, boundary);
        case 'none'
            % assume that data are already unwrapped 
        otherwise
            error('Missing method')
    end
end
