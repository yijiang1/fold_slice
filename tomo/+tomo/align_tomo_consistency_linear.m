%    ALIGN_TOMO_CONSISTENCY_LINEAR   self consitency based alignment procedure based on the ASTRA toolbox 
%       it can align both in horizontal and vertical dimension 
%
%   [optimal_shift,par, rec, err] = ...
%             align_tomo_consistency_linear(sinogram_0,weights, angles, Npix, optimal_shift, par, varargin )
%
% Inputs: 
%   **sinogram_0        - real value sinogram (ie phase difference or unwrapped phase), see "unwrap_data_method" for more details 
%   **angles            - angle in degress
%   **Npix              - size of the reconstructed field 
%   **optimal_shift     - initial guess of the shift 
%   **par               - parameter structure 
%  *optional*
%   **align_vertical	- (bool) allow vertical alignment 
%   **align_horizontal  - (bool) allow horizontal alignment 
%   **mask_threshold    - (scalar) values < threshold will be considered empty, [] == auto guess 
%   **apply_positivity  - (bool) remove negative values 
%   **high_pass_filter  - (scalar) bandpass filter applied on the data, 1 = no filtering, eps = maximal filtering, common value ~0.01 to avoid low spatial freq. errors 
%   **max_iter - (int) maximal number of iterations for the alignment 
%   **lamino_angle - (scalar) tilt angle of the tomographic axis 
%   **show_projs - (bool) show a movie of reprojections and input
%   OTHER INPUTS AND THEIR DEFAULT VALUES ARE DESCRIBED IN CODE
% 
% *returns* 
%   ++optimal_shift         - optimal shift of the inputs in order to maximize the tomographic consistency 
%   ++par                   - tomographic parameter structure, keeps updated values of the geometric parameters 
%   ++rec                   - final reconstruction given the selected binning
%   ++err                   - evolution of misfit between the measured and the modelled projections 


function [optimal_shift,params, rec, err] = ...
            align_tomo_consistency_linear(sinogram,weights, angles, Npix, optimal_shift, params, varargin )
    


    import tomo.*
    import utils.*
    import math.*
    utils.verbose(struct('prefix', 'align'))


    parser = inputParser;
    % general  parameters 
    parser.addParameter('align_vertical',  true , @islogical )
    parser.addParameter('align_horizontal', false , @isnumeric )
    parser.addParameter('high_pass_filter',  0.02 , @isnumeric ) % bandpass filter applied on the data, 1 = no filtering, eps = maximal filtering, common value ~0.01
    parser.addParameter('min_step_size',  0.01 , @isnumeric )  % minimum step size, optimization is ended once step is lower 
    parser.addParameter('max_iter',  50 , @isint )          % maximal number of iterations executed if min_step_size is not reached
    parser.addParameter('verbose',  1 , @isnumeric ) % change verbosity of the code 
    parser.addParameter('is_laminography',  false , @isnumeric ) % assume that reconstruction geometry is laminography
    parser.addParameter('is_interior_tomo',  false , @isnumeric ) % assume that reconstruction geometry is interior tomography

    parser.addParameter('refine_geometry', false, @islogical ) % try to refine global errors in geometry - lamino angle, tilt, skewness
    parser.addParameter('refine_geometry_parameters', {'shear_angle', 'tilt_angle', 'lamino_angle','asymmetry'}, @iscell ) % try to refine global errors in geometry - lamino angle, tilt, skewness
    parser.addParameter('use_GPU',  gpuDeviceCount > 0 , @islogical) % apply affine deformation on the projections 
    parser.addParameter('momentum_acceleration',  false , @islogical) % accelerate convergence by momentum gradient descent method 
    parser.addParameter('online_tomo',  false , @islogical)  % if true, dont expect any feedback from user  

    % extra contraints 
    parser.addParameter('mask_threshold',  [] , @isnumeric ) % threshold for mask estimation 
    parser.addParameter('use_mask',  false , @islogical )    % autoestimate mask and apply 
    parser.addParameter('use_localTV',  false , @islogical )  % additional regularization
    parser.addParameter('localTV_lambda',  1e-4 , @isnumeric )  % added by YJ. parameter so control TV strength
    parser.addParameter('apply_positivity',  false , @islogical )
    parser.addParameter('apply_soft_threshold',  false , @islogical ) % added by YJ
    parser.addParameter('soft_threshold',  1e-4 , @isnumeric )  % added by YJ

    % geometry parameters 
    parser.addParameter('lamino_angle',  90 , @isnumeric )   % laminography title (with respect to the beam ))
    parser.addParameter('tilt_angle',  0 , @isnumeric )     % rotation of camera around axis of the beam 
    parser.addParameter('skewness_angle',  0 , @isnumeric ) % skewness of the projection coordinate axis 
    parser.addParameter('pixel_scale',  [1,1] , @isnumeric ) % ratio between the vertical and horizontal axis 
    parser.addParameter('CoR_offset', [] , @isnumeric ) % offset of the center of rotation, default is center of projection 
    parser.addParameter('CoR_offset_v', [] , @isnumeric ) % offset of the center of rotation, default is center of projection 

    % reconstruction method parameters 
    parser.addParameter('deformation_fields', [])  % assume deformated sample and use these fields
    parser.addParameter('inv_deformation_fields', [])  % get also inversion, otherwise use -deformation_fields    parser.addParameter('plot_results',  true , @islogical ) % plot results 
    parser.addParameter('filter_type',  'ram-lak' , @isstr ) % change verbosity of the code 
    parser.addParameter('freq_scale',  '1' , @isnumeric ) % change verbosity of the code 
    parser.addParameter('plot_results',  true , @islogical ) % plot results 
    parser.addParameter('show_projs',  false , @islogical ) % show a movie of reprojections and input 

    % data related parameters 
    parser.addParameter('binning',  1 , @isint )            % downsample dataset by binning to make alignment faster and more robust 
    parser.addParameter('unwrap_data_method',  'fft_1D' , @isstr )  % options: "none" (inputs is directly phase or amplitude), "fft_1d" or "diff" assume  that inputs is phase derivative 
    parser.addParameter('valid_angles', [], @islogical ) % bool array of valid projections used for refienement  
    parser.addParameter('selected_roi',  {} , @iscell ) % field of view considered for alignment for laminography
    parser.addParameter('vert_range',  [] , @isnumeric) % vertical range considered for alignment for standard tomo 
    parser.addParameter('affine_matrix',  [] , @isnumeric) % apply affine deformation on the projections 

    
    % plotting 
    parser.addParameter('plot_results_every',  10 , @isnumeric) % plot results every N seconds 
    parser.addParameter('position_update_smoothing',  0 , @isnumeric) % enforce smoothness of the updates, useful in the initial alignment 
    parser.addParameter('windowautopos',  true , @isnumeric) % automatically place the plotted windows 

    parser.parse(varargin{:})
    r = parser.Results;

    % load all to the param structure 
    par = params; 
    for name = fieldnames(r)'
        if ~isfield(par, name{1}) || ~ismember(name, parser.UsingDefaults) % prefer values in param structure 
            par.(name{1}) = r.(name{1});
        end
    end
    
    Nangles = length(angles); 
    assert(isa(sinogram, 'single') || isa(sinogram, 'uint16'), 'Input sinogram has to be single precision')
    assert(size(sinogram,3) == Nangles, 'Number of angles does not correspond to number of projections')
    assert(size(optimal_shift,1) == Nangles, 'Number of angles does not correspond to dimensions of optimal_shift')
    assert(length(par.valid_angles) == Nangles, 'Number of angles does not correspond to dimensions of par.valid_angles')
    if par.refine_geometry && ~(par.align_horizontal || par.align_vertical)
       warning('if par.refine_geometry = true, it is recommended to set par.align_horizontal =  true, par.align_vertical = true') 
    end
    verbose(1,'Starting align_tomo_consistency_linear')
    [Nlayers,Nw,Nangles] = size(sinogram);
    tic
    warning on 
    
    if abs(mean(par.lamino_angle) - 90) > 1 && ~par.is_laminography
        warning('Use par.is_laminography == true for lamino_angle < 90')
    end
    if isempty(par.valid_angles)
        valid_angles = true(Nangles,1); 
    else
        valid_angles = par.valid_angles; 
    end
    
    if ~par.is_laminography && isempty(par.selected_roi )
        % find optimal vertical range 
        if isempty(par.vert_range)
            par.vert_range = [1, Nlayers]; 
        end
        vrange0 =  [max(par.vert_range(1),round(1+max(optimal_shift(:,2)))),min(par.vert_range(end), floor(Nlayers + min(optimal_shift(:,2))))]; 
        vrange(2) = min(vrange0(2),Nlayers);
        vrange(1) = max(1,vrange0(1));
        
        % help with splitting in ASTRA (GPU memory limit)
        vrange_center = ceil(mean(vrange)); 
        estimated_split_factor = 2^nextpow2(ceil(prod(Npix(1)*Npix(min(2,end))*Nlayers)*8/1e9/par.binning^3)*par.binning); 
        Nvert = floor((vrange(2)-vrange(1)+1) / estimated_split_factor)*estimated_split_factor;
        vrange(1) = ceil(vrange_center - Nvert/2); 
        vrange(2) = floor(vrange_center + Nvert/2-1); 


        if any(vrange ~= vrange0)
           verbose(1,'Changing vertical range from %i:%i to %i:%i', vrange0(1),vrange0(2), vrange(1), vrange(2)) 
        end
        if isempty(vrange(1):vrange(2)) || (par.align_vertical && length(vrange(1):vrange(2)) < 10*par.binning)
            error('Too small par.vert_range, extend the alignment range'); 
        end 
        %% limit the vertical range
        ROI = {vrange(1):vrange(2), ':'}; 
        Nlayers = length(ROI{1});
    else
        ROI = par.selected_roi; 
    end
    
    if all(weights(:)==1)
        weights = [];
    end
    
    if isempty(par.affine_matrix)
       par.affine_matrix = diag([1,1]); 
    end
   
    if ~par.align_horizontal && ~par.align_vertical
        warning('align_horizontal and align_vertical are both false')
    end
    
    if ~isempty(par.deformation_fields)
        % interpolate the deformation fields to fit the reduced vertical range 
        for jj = 1:numel(par.deformation_fields)
            Nl_deform = size(r.deformation_fields{jj}{1},3); 
            vrange_tmp=0.5+vrange / Nlayers*Nl_deform;
            % oversample the field twice in direction of the cropped axis 
            interp_ax = min(Nl_deform,max(1,linspace(vrange_tmp(1),vrange_tmp(2), 2*(ceil(vrange_tmp(2))- floor(vrange_tmp(1)))))); 
            for ii = 1:numel(par.deformation_fields{jj})
                par.deformation_fields{jj}{ii} = permute(interp1(permute(par.deformation_fields{jj}{ii},[3,1,2]),interp_ax),[2,3,1]);
                par.inv_deformation_fields{jj}{ii} = permute(interp1(permute(par.inv_deformation_fields{jj}{ii},[3,1,2]),interp_ax),[2,3,1]);
                par.deformation_fields{jj}{ii} = par.deformation_fields{jj}{ii}/par.binning;
                par.inv_deformation_fields{jj}{ii} = par.inv_deformation_fields{jj}{ii}/par.binning;
            end
        end
    end

    verbose(1,'Shifting sinograms and binning = %i', par.binning) % + crop them into smaller field of view

    % depending on the unwrapping method, te fft interpolation needs different shift to account for
    % "binning/downsampling" of the data 
    if strcmpi(par.unwrap_data_method, 'none')
        interp_sign = -1 ;
    else
        interp_sign = 1 ;
    end
    
    % define parameters structure  that determine behaviour of the
    % tomo.block_fun function 
    Np_sinogram = size(sinogram); 
    param_struct = struct('GPU_list',par.GPU_list, 'full_block_size', Np_sinogram, 'use_fp16', false); 

    % shift to the last optimal position  + remove edge issues 
    sinogram =  tomo.block_fun(@imshift_generic, sinogram, optimal_shift,Np_sinogram,par.affine_matrix,5,ROI,par.binning, 'fft',interp_sign, param_struct);
    % shift also the  weights to correspond to the sinogram
    if ~isempty(weights)
        if isa(weights, 'uint8')
            weights = single(weights)/255; % load the weights from uint8 format
        end
        assert(all(sum(sum(weights)) > 0), sprintf('Provided "weights" contain %i projections with empty mask', sum(sum(sum(weights))==0)))
        weights =  tomo.block_fun(@imshift_generic, single(weights), optimal_shift,Np_sinogram,par.affine_matrix,0, ROI,par.binning, 'linear', param_struct);
    end

    % sort by angle if requested, but only after binning to make it faster
    if par.showsorted
        [~,ang_order] = sort(angles);
        % sort / crop the angles , slighly improves speed of ASTRA 
        sinogram = sinogram(:,:,ang_order); 
        try weights = weights(:,:,ang_order); end
        angles = angles(ang_order); 
        optimal_shift = optimal_shift(ang_order,:); 
        valid_angles = valid_angles(ang_order);
    else
        ang_order = 1:Nangles;
    end
    par.ang_order = ang_order; 
    par.inv_ang_order(ang_order) = 1:Nangles;
    affine_matrix = [];
 
    %% %%%%%%%%%%%%%%%% initialize astra %%%%%%%%%%%%%%%%
    [Nlayers,width_sinogram,~] = size(sinogram);
    if gpuDeviceCount
        %% %%%%%%%%%%  initialize GPU %%%%%%%%%%%%%%%
        gpu  = gpuDevice();
        if ~isempty(par.GPU_list) && gpu.Index ~= par.GPU_list(1)
            gpu  = gpuDevice(par.GPU_list(1));
        end

        % keep small datasets fully in GPU to speed it up 
        keep_on_GPU =  (gpu.AvailableMemory >  numel(sinogram) * 4 * 8 + Npix(1)*Npix(min(2,end))*Nlayers*4*8/par.binning^3) && ...
                         (numel(sinogram) < intmax('int32')) && ...
                         (prod(Npix) < intmax('int32')); 

        % split it among multiple GPU only for larger datasets (>4GB)
        if ~(numel(sinogram) * 4 > 1e9 && length(par.GPU_list) > 1  && ~keep_on_GPU)
            par.GPU_list = gpu.Index;
        end

        block_cfg = struct('GPU_list',par.GPU_list); 
    else
        block_cfg = struct();
        keep_on_GPU = false; 
    end
    %keep_on_GPU = false; 
    block_cfg.verbose_level = utils.verbose();
    
    if  ~keep_on_GPU
        reset(gpuDevice)  % better reset the GPU before processing large datasets
    end
    
    if ~par.use_GPU
        par.tomo_solver = @FBP_CPU ; 
        par.padding = 0; 
    elseif ~par.is_laminography && length(par.GPU_list) <= 1
        par.tomo_solver = @FBP_zsplit ; 
        par.padding = 0; % padding method in when the filter is applied 
    else
        par.tomo_solver = @FBP; 
        par.padding = 'symmetric'; % padding method in when the filter is applied 
    end

    %% prepare some auxiliary variables 
    win = tukeywin(width_sinogram, 0.2)';  % avoid edge issues 
    if Nlayers > 10 && (par.align_vertical ); win = tukeywin(Nlayers, 0.2)*win; end 
        
    shift_upd_all = nan(par.max_iter,Nangles,2, 'single');
    shift_upd_all(1,:,:) = 0;
    shift_velocity =  zeros(Nangles,2, 'single');
    shift_total =  zeros(Nangles,2, 'single');
    err = nan(1,Nangles,'single'); 
    
    time_plot = tic; 
    last_round = false; 

    % ASTRA needs the reconstruction to be dividable by 32 othewise there
    % will be artefacts in left corner 
    Npix = ceil(Npix/par.binning);
    if isscalar(Npix)
        Npix = [Npix, Npix, Nlayers];
    elseif length(Npix) == 2
        Npix = [Npix, Nlayers];
    end
    
    % !! important for binning => account for additional shift of the center
    % of rotation after binning, for binning == 1 the correction is zero
    par.rotation_center = [Nlayers, width_sinogram]/2; 

    if ~isempty(par.CoR_offset)
        par.rotation_center(2) = par.rotation_center(2) + par.CoR_offset/par.binning; 
    end
    if ~isempty(par.CoR_offset_v)
        par.rotation_center(1) = par.rotation_center(1) + par.CoR_offset_v/par.binning; 
    end
    
    %%  weights of the pixels in the sinograms, important for
    %%  laminography or interior tomography
    if isempty(weights); weights = 1; end
    %[~,circulo]=utils.apply_3D_apodization(zeros(Npix),0,0, 5); %old code
    apodize = 0;
    if isfield(par,'apodize')
        apodize = par.apodize;
    end
    circulo = utils.get_apodization_mask(Npix, apodize/par.binning);   %modified by YJ for tomography of flat object
    %circulo = 1; %modified by YJ to prevent error
    %figure
    %imagesc(circulo)
    % if the array is small enough, move the calculations fully on GPU 
    if keep_on_GPU 
        %keep all on gpu only for small arrays 
        sinogram = Garray(sinogram);
        weights  = Garray(weights);
        circulo  = Garray(circulo);
    end
    
    %% geometry parameters structure 
    geom.tilt_angle = 0; 
    geom.skewness_angle = 0; 
    geom.asymmetry = 0; 
    geom.lamino_angle = par.lamino_angle;

    time_per_iteration = nan; 
    for ii = 1:par.max_iter
        t0 = tic; 
        if verbose() > 0
            verbose(1,'Iteration %i / %i\ttime: %4.2gs\tError %7.5s', ii, par.max_iter, time_per_iteration, mean(err(max(1,ii-1),:)))
        else
            utils.progressbar(ii, par.max_iter) 
        end
        sinogram_shifted = sinogram; 
        weights_shifted   =  weights; 

        %% shift the sinogram and weights 
        affine_matrix  =    compose_affine_matrix(1, geom.asymmetry, -geom.tilt_angle, -geom.skewness_angle);
        par.lamino_angle =  geom.lamino_angle;
        sinogram_shifted =  tomo.block_fun(@utils.imdeform_affine_fft, sinogram_shifted, affine_matrix, shift_total, block_cfg);

        if ~ismatrix(weights_shifted)
            % use linear shift in order to prevent periodic boundary issues 
            % !! still very slow operation
            weights_shifted =   tomo.block_fun(@imshift_linear, weights_shifted, shift_total(:,1),shift_total(:,2), 'linear', block_cfg);
        end

        weights_shifted = max(0,bsxfun(@times, weights_shifted, win)) ; % apply filter to avoid edge issues 

        
        if ~strcmpi( par.unwrap_data_method, 'none')
            sinogram_shifted = tomo.block_fun(@unwrap_data, sinogram_shifted,  par.unwrap_data_method, par.air_gap/par.binning, block_cfg); 
        end
        if ii == 1; MASS = gather(median(tomo.block_fun(@(x)mean2(abs(x)), sinogram_shifted, struct('use_GPU', false, 'verbose_level', utils.verbose())))); end

        %% %%%%%%%%%%%%  FBP method  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        verbose(2,'FBP')
        try
            [rec,cfg,vectors] = get_reconstruction(sinogram_shifted, angles, Npix, keep_on_GPU, par);
        catch myErr
            if length(par.GPU_list)>1, delete(gcp('nocreate')); end
            disp(getReport(myErr))
            if ~par.online_tomo
                keyboard
            else
               rethrow(myErr) 
            end
        end
        %% %%%%%%%%% apply all available prior knowledge %%%%%%%%%%%%%%%%%%%%%%%
        if par.use_mask
            if mod(ii,5) == 1 
               mask = get_mask(rec, par, circulo); 
            end
            % apply mask 
            rec = rec .* mask;
        else
            rec = rec .* circulo; % apply at least always mask on the out of view regions
                                  % it seems that stability of the method is greatly improved if at
                                  % least very loose mask is applied 
        end
        
        if par.use_localTV
            rec = regularization.local_TV3D_chambolle(rec, par.localTV_lambda, 10);
            %rec = regularization.local_TV3D_grad(rec, par.localTV_lambda, 10);
        end
        if par.apply_positivity
            %rec = max(0,rec);  % significantly improves the alignment stability , dont use for laminography / missing wedge tomo 
            %rec = max(par.soft_threshold,rec);  % significantly improves the alignment stability , dont use for laminography / missing wedge tomo 
            rec(rec<par.soft_threshold) = 0;
        end
        if par.apply_soft_threshold
           rec = sign(rec).*max(abs(rec)-par.soft_threshold,0); 
        end
        %% %%%%%%%%%%%% center reconstruction or keep initial position %%%%%%%%%
        % important to remove 
        if par.align_horizontal && ~par.is_laminography && ~par.is_interior_tomo  % do not center it in the first step
            [x,y,mass] = center(sqrt(max(0,rec))+eps); % abs seems to be more stable than max(0,x) even for missing wedge or laminography 
            if Npix(3) > 5; mass = mass .* reshape(tukeywin(Npix(3),0.1),1,1,[]); end
            % more robust estimation of center 
            rec_center(1) = gather(mean(x.*mass)./mean(mass)); 
            rec_center(2) = gather(mean(y.*mass)./mean(mass)); 

            if ii == 1 
                if par.center_reconstruction
                    rec_center_0 = [0,0] ;  % ideal rotation center
                else
                    rec_center_0 = rec_center;  % just keep the position of the first iteration 
                end
                verbose(2,'Reconstruction center %g %g \n', rec_center_0)
            end
            
            % shift direction, go slowly
            shift_rec = -0.5*(rec_center - rec_center_0); 
             
            if par.center_reconstruction
                verbose(2,'Centering reconstruction');
            else
                verbose(2,'Keeping center of mass of reconstruction');
            end

            % avoid drifts of the reconstructed volume 
            rec = tomo.block_fun(@imshift_fft, rec, shift_rec(1), shift_rec(2), block_cfg); 
            if par.use_mask
                mask = tomo.block_fun(@imshift_fft, mask, shift_rec(1), shift_rec(2), block_cfg) > 0.5; 
            end
        end

        %% %%%%%%%%  Get computed sinogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        verbose(2,'Forward projection')
        try
            sinogram_model = get_projections(rec,cfg, vectors,par);    
        catch myErr
            if length(par.GPU_list)>1, delete(gcp('nocreate')); end
            disp(getReport(myErr))
            if ~par.online_tomo
                keyboard
            else
               rethrow(myErr) 
            end
        end

        
        %% find shift of the data sinogram to match the "optimal sinogram"
        verbose(2,'Find optimal shift');
        if keep_on_GPU
            %% if possible, move to GPU 
            sinogram_model = Garray(sinogram_model); 
            sinogram_shifted = Garray(sinogram_shifted); 
        else
            sinogram_model = gather(sinogram_model); 
        end
        
        %% refine geometry 
        if par.refine_geometry && gpuDeviceCount > 0 
            [geom] = refine_geometry(rec, sinogram_shifted,sinogram_model,weights_shifted, angles, ii,Npix, Nlayers,width_sinogram, Nangles,par,geom);
        end

        %% FOR DEBUGGING - play movie of the synthetic and real projection after filtering 
        if par.show_projs
            plotting.smart_figure(1234)
            plotting.imagesc3D(cat(1, sinogram_model .* weights_shifted, sinogram_shifted .* weights_shifted)); 
            axis off image ; colormap bone; grid on
            title('Top: model  Bottom:  filtered data')
        end
        %% find optimal shift 
        try
            [shift_upd, err(ii,:)] = tomo.block_fun(@find_optimal_shift,sinogram_model, sinogram_shifted,weights_shifted, MASS, par, block_cfg); 
        catch err
            disp(getReport(err))
            keyboard
        end
        clear sinogram_model
        
        % do not allow more than 0.5px per iteration !!
        shift_upd = min(0.5, abs(shift_upd)).*sign(shift_upd)*par.step_relaxation;

        % store update history for momentum gradient acceleration  
        shift_upd_all(ii,:,:) = reshape(shift_upd, [1,Nangles,2]);
                
        %% find correlation between the subsequent updates -> accelerate convergence, something like momentum method 
        % -> speed up convergence
        if par.momentum_acceleration
            momentum_memory = 2; 
            max_update = quantile(abs(shift_upd(valid_angles,:)), 0.995); 
            if ii > momentum_memory
                [shift_upd, shift_velocity]= add_momentum(shift_upd_all(ii-momentum_memory:ii,:,:), shift_velocity, max_update*par.binning < 0.5 );
                %shift_upd_all(ii,:,:) = reshape(shift_upd, [1,Nangles,2]);
            end
        else
            if ii > 1
                verbose(1, 'Correlation between two last updates x:%.3f%% y:%.3f%%', 100*corr(squeeze(shift_upd_all(ii-1,:,1))', squeeze(shift_upd_all(ii,:,1))'), 100*corr(squeeze(shift_upd_all(ii-1,:,2))', squeeze(shift_upd_all(ii,:,2))'))
            end
        end


        shift_upd(:,2) = shift_upd(:,2) - median(shift_upd(:,2));
        
        % prevent outliers when the code decides to quickly oscilate around the solution
        max_step = min(quantile(abs(shift_upd), 0.99), 0.5); 
        % do not allow more than 0.5px per iteration (multiplied by binning factor) !!
        shift_upd = min(max_step, abs(shift_upd)).*sign(shift_upd); 
        
            
        % remove degree of freedom in the vertical dimension (avoid drifts)     
         if par.align_horizontal && par.is_laminography
             orthbase = [sind(angles(:)), cosd(angles(:))]; % 
             coefs = (orthbase'*orthbase) \ (orthbase'*shift_upd(:,1));
             % avoid object drifts within the reconstructed FOV
             rigid_shift = orthbase*coefs; 
             shift_upd(:,1) = shift_upd(:,1) - rigid_shift(:,1);
         end   

        % update the total position shift 
        shift_total = shift_total + shift_upd;

            
        % remove degree of freedom in the vertialignment_ROIcal dimension (avoid drifts)     
%          if ~( par.center_reconstruction) && par.is_laminography
%              orthbase = [sind(angles(:)), cosd(angles(:)) ,ones(Nangles,1)]; % 
%              coefs = (orthbase'*orthbase) \ (orthbase'*shift_total);
%              % avoid object drifts within the reconstructed FOV
%              rigid_shift = orthbase*coefs; 
%              shift_total(:,1) = shift_total(:,1) - rigid_shift(:,1);
%          end   

        % enforce smoothness of the estimate position update -> in each
        % iteration smooth the accumulated position update 
        % this helps againts discontinuities in the update 
        if par.position_update_smoothing
            for kk = 1:2
                shift_total(:,kk) =  smooth(shift_total(:,kk), max(0, min(1, par.position_update_smoothing)) * Nangles); 
            end
        end
        
        max_update = max(quantile(abs(shift_upd(valid_angles,:)), 0.995)); 
        verbose(1,'Maximal step update: %4.2g px stopping criterion: %4.2g px ', ...
            max_update * par.binning, par.min_step_size );
        if max_update * par.binning < par.min_step_size        
            verbose(1,'Minimal step limit reached')
            last_round = true;  % stop iterating if converged 
        end

        if par.plot_results  &&  (toc(time_plot) > par.plot_results_every || last_round || ii == par.max_iter ) % avoid plotting in every iteration, it is too slow ... 
        %%%%%%%%%%%%%%%%%%%% plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            plot_alignment(rec, sinogram_shifted, weights_shifted, err,shift_upd,shift_total, angles, valid_angles, ii, par)
            time_plot = tic; 
        end
        
        if last_round
            utils.progressbar(par.max_iter, par.max_iter) 
            break
        end
        time_per_iteration = toc(t0); 
        clear sinogram_shifted
    end

    if ~par.is_interior_tomo
        % vertical offset is a degree of freedom => minimize of offset 
        shift_total(:,2) = shift_total(:,2) - median(shift_total(:,2));  
    end

    %% prepare outputs to be exported     
    % sort / crop the angles 
    optimal_shift(par.ang_order,:) = optimal_shift + (shift_total )* par.binning;
    err(:,par.ang_order) = err; 
    
    verbose(1,'Done')
    params.center_reconstruction = false;
    params.lamino_angle          = geom.lamino_angle; 
    params.tilt_angle            = par.tilt_angle     +geom.tilt_angle; 
    params.skewness_angle        = par.skewness_angle +geom.skewness_angle; 
    params.pixel_scale           = par.pixel_scale; 

    rec = gather(rec);

    
    utils.verbose(struct('prefix', 'template'))

end


function geom = refine_geometry(rec, sinogram_measured,sino_model,weights, angles , iter, Npix, Nlayers,width_sinogram, Nangles,par,geom)
    import math.*
    import utils.*
    persistent geometry_corr
    if iter == 1
        geometry_corr = [mean(par.lamino_angle), mean(geom.tilt_angle), mean(geom.skewness_angle),mean(geom.asymmetry,1)]; 
    end 
    
    gpu = gpuDevice; 
    keep_on_GPU =  (gpu.AvailableMemory >  numel(sino_model) * 4 * 16) && ...
                    (numel(sino_model) < intmax('int32')) && ...
                    (prod(Npix) < intmax('int32')); 
                
    if keep_on_GPU
        sino_model = Garray(sino_model); 
        sinogram_measured = Garray(sinogram_measured); 
        GPU_list = par.GPU_list(1);  % use only single GPU if keep_on_GPU is used 
    else
        GPU_list = par.GPU_list;
    end
    
    resid_sino = get_resid_sino(sino_model, sinogram_measured, par.high_pass_filter);

    step_relaxation = 0.01; 
    [dX,dY] = math.get_img_grad(sino_model); 
    lamino_angle_plot = [];
    tilt_angle_plot = [];
    skewness_angle_plot = [];
    if par.is_laminography && any(ismember(par.refine_geometry_parameters, 'lamino_angle'))
        % get laminography angle correction 
        optimal_shift = get_geometry_corr(rec,resid_sino,weights, Npix, Nlayers, width_sinogram, angles,[1,0,0,0],GPU_list, par);
        lamino_angle_plot = geom.lamino_angle +  step_relaxation*gather((optimal_shift));
        geom.lamino_angle = geom.lamino_angle +  step_relaxation*gather(median(optimal_shift)); 
    end
        
    % get tilt angle correction 
    if any(ismember(par.refine_geometry_parameters, 'tilt_angle'))
        Dvec = (dX .* linspace(-1,1,size(dX,1))' - dY .* linspace(-1,1,size(dY,2))); 
        optimal_shift = get_GD_update(Dvec, resid_sino, weights,  par.high_pass_filter); 
        tilt_angle_plot = geom.tilt_angle +step_relaxation*rad2deg(gather(optimal_shift));
        geom.tilt_angle = geom.tilt_angle +step_relaxation*rad2deg(gather(optimal_shift));
    end
    
    % get shear gradient
    if any(ismember(par.refine_geometry_parameters, 'shear_angle'))
        Dvec = dY .* linspace(-1,1,size(dY,2)); 
        optimal_shift = get_GD_update(Dvec, resid_sino, weights,  par.high_pass_filter); 
        skewness_angle_plot = geom.skewness_angle +step_relaxation*rad2deg(gather(optimal_shift));

        geom.skewness_angle = geom.skewness_angle +step_relaxation*rad2deg(gather(optimal_shift));
    end
    
    % get pixels err
    if any(ismember(par.refine_geometry_parameters, 'asymmetry'))
        Dvec = dX .* linspace(-1,1,size(dX,2)); 
        optimal_shift_1 = get_GD_update(Dvec, resid_sino, weights,  par.high_pass_filter); 
        Dvec = dY .* linspace(-1,1,size(dX,1))'; 
        optimal_shift_2 = get_GD_update(Dvec, resid_sino, weights,  par.high_pass_filter); 
        optimal_shift = [optimal_shift_1, optimal_shift_2]; 
        if ~par.is_laminography; optimal_shift = optimal_shift - mean(optimal_shift); end % remove extra degrees of freedom
        geom.asymmetry = gather(optimal_shift_2./optimal_shift_1); % ignote difference in pixel scale
    end

    geometry_corr(end+1,:) = [mean(geom.lamino_angle), mean(geom.tilt_angle), mean(geom.skewness_angle), mean(geom.asymmetry)]; 

    if par.is_laminography
         geom.lamino_angle = mean(geom.lamino_angle); 
         %modified by YJ
         geom.skewness_angle = mean(geom.skewness_angle); 
         %geom.skewness_angle = (geom.skewness_angle); 
         geom.tilt_angle = mean(geom.tilt_angle); 
         %geom.tilt_angle = (geom.tilt_angle); 
    end
    
    if mod(iter, 5) == 1
        plotting.smart_figure(123)
        subplot(2,4,1)
        plot(geometry_corr(:,1));
        set(gca, 'xscale', 'log')
        ylabel('AVG Laminography angle');
        grid on; axis tight 
        xlabel('Iteration')
        subplot(2,4,2)
        plot( mean(par.tilt_angle) + geometry_corr(:,2));
        set(gca, 'xscale', 'log')
        ylabel('AVG Projection rotation');
        grid on; axis tight 
        xlabel('Iteration')
        subplot(2,4,3)
        plot( mean(par.skewness_angle) + geometry_corr(:,3));
        set(gca, 'xscale', 'log')
        ylabel('AVG Projection skewness') ; 
        grid on; axis tight 
        xlabel('Iteration')
        subplot(2,4,4)
        plot(geometry_corr(:,4));
        set(gca, 'xscale', 'log')
        ylabel('AVG Relative pixel scale') ; 
        grid on; axis tight 
        xlabel('Iteration')

        subplot(2,4,5)
        plot(lamino_angle_plot);
        ylabel('Laminography angle');
        grid on; axis tight 
        xlabel('Slice')
        subplot(2,4,6)
        plot(par.tilt_angle+tilt_angle_plot);
        ylabel('Projection rotation');
        grid on; axis tight 
        xlabel('Slice')
        subplot(2,4,7)
        plot(par.skewness_angle+skewness_angle_plot );
        ylabel('Projection skewness') ; 
        grid on; axis tight 
        xlabel('Slice')
        subplot(2,4,8)
        plot(geom.asymmetry);
        ylabel('Relative pixel scale') ; 
        grid on; axis tight 
        xlabel('Slice')
        plotting.suptitle('Evolution of the geometry refinement')
        drawnow
    
    end
    
end

function optimal_shift = get_GD_update(dX, resid, weights, filter)
    % auxiliary function that calculate the optimal step length for optical
    % flow based image alignment 
    import math.* 
    dX = imfilter_high_pass_1d(dX, 2,filter);     
    optimal_shift = squeeze(sum2(weights .* resid .* dX) ./ sum2(weights .* dX.^2));
end
        
function optimal_shift = get_geometry_corr(rec, resid_sino,weights,Npix, Nlayers, width_sinogram, angles,search_dim, GPU_list,par)
    % gradient descent base calculation of optimal geometry correction 
    
    import math.* 
    par.GPU_list = GPU_list; 
    delta = 0.01; % small step used to calculate numerically the gradient 
    [cfg, vectors] = ...
        astra.ASTRA_initialize(Npix, [Nlayers, width_sinogram],angles,par.lamino_angle-delta*search_dim(1),par.tilt_angle-delta*search_dim(2), par.pixel_scale-[delta*search_dim(4),0], par.rotation_center, delta*search_dim(3)); 
    sinogram{1} = get_projections(rec,cfg, vectors,par); 
    [cfg, vectors] = ...
        astra.ASTRA_initialize(Npix, [Nlayers, width_sinogram],angles,par.lamino_angle+delta*search_dim(1),par.tilt_angle+delta*search_dim(2), par.pixel_scale+[delta*search_dim(4),0], par.rotation_center, delta*search_dim(3)); 
    sinogram{2} = get_projections(rec,cfg, vectors,par); 
    dsino = (sinogram{1}-sinogram{2})/(2*delta); 
    dsino = imfilter_high_pass_1d(dsino, 2, par.high_pass_filter);
    optimal_shift = squeeze(sum2(weights .* resid_sino .* dsino) ./ sum2(weights .* dsino.^2));

end

function [rec, cfg, vectors] = get_reconstruction(sinogram, angles, Npix, keep_on_GPU, par)
    % auxiliar function to get FBP tomography reconstruction given the provided parameters 
    % automatically apply evolving deformation field if provided 

    import tomo.*
    import utils.verbose

    [Nlayers, width_sinogram,Nangles] = size(sinogram);
    N_GPU = max(1,length(par.GPU_list)); 

    [cfg, vectors] = ...
        astra.ASTRA_initialize(Npix, [Nlayers, width_sinogram],angles,par.lamino_angle,par.tilt_angle, par.pixel_scale, par.rotation_center, par.skewness_angle); 

    % find optimal split of the dataset for given GPU 
    split = astra.ASTRA_find_optimal_split(cfg, N_GPU,1, 'back');

    params = {'valid_angles',par.valid_angles, ...
            'GPU', par.GPU_list, 'verbose', par.verbose_level > 2, 'keep_on_GPU', keep_on_GPU, ...
            'filter', par.filter_type, 'filter_value', par.freq_scale, ...
            'use_derivative', strcmpi( par.unwrap_data_method, 'diff'), 'padding', par.padding}; 

    if isempty(par.inv_deformation_fields)
        rec  = par.tomo_solver(sinogram, cfg, vectors, split, params{:});
    else
        rec = 0; 
        Nblocks = length(par.inv_deformation_fields); 
        Bsize = ceil(Nangles / Nblocks); 
        for ll = 1:Nblocks
             ids = par.inv_ang_order(1+(ll-1)*Bsize:min(Nangles, ll*Bsize)); 
             rec  = rec +1/Nblocks* FBP(sinogram, cfg, vectors, split,...
                     'deformation_fields', par.inv_deformation_fields{ll}, params{:}, ...
                     'valid_angle', ids);
        end
    end
end

function sinogram = get_projections(rec,cfg, vectors, par)
    % auxiliar function to get tomography projections given the provided parameters 
    % automatically apply evolving deformation field if provided 
    
    import utils.verbose
    
    N_GPU = max(1,length(par.GPU_list)); 
    cfg.iProjAngles = size(vectors,1); 
    split = astra.ASTRA_find_optimal_split(cfg, N_GPU,1, 'fwd');
    
    if ~par.use_GPU
         sinogram = tomo.radon_wrapper(rec,cfg, vectors); 
    elseif par.is_laminography && N_GPU == 1
        % single GPU code, seems to be faster for laminography 
        sinogram = astra.Ax_partial(rec,cfg, vectors,split,...
             'GPU', par.GPU_list,'verbose', verbose());
    elseif isempty(par.inv_deformation_fields)
        % multiGPU parallel code, for normal tomo as fast as Ax_partial
        % but more memory efficient 
        sinogram = tomo.Ax_sup_partial(rec,cfg, vectors,[1, 1,N_GPU] ,...
            'GPU', par.GPU_list, 'split_sub', split, 'verbose' , verbose());
    else
        % deformation tomography, solve blockwise 
        Nblocks = length(par.inv_deformation_fields); 
        Bsize = ceil(cfg.iProjAngles / Nblocks);
        sinogram = zeros(cfg.iProjV,cfg.iProjU,cfg.iProjAngles, 'like', rec); 
        for ll = 1:Nblocks
            ids = par.inv_ang_order(1+(ll-1)*Bsize:min(cfg.iProjAngles, ll*Bsize)); 
            sinogram(:,:,ids) = tomo.Ax_sup_partial(rec,cfg, vectors(ids,:),[1,1,N_GPU] ,...
            'GPU', par.GPU_list, 'split_sub', split, 'verbose', verbose(),  'deformation_fields', par.deformation_fields{ll});
        end
    end
end

function [shift,velocity_map] = add_momentum(shifts_memory, velocity_map,acc_axes)
    % function for accelerated momentum gradient descent. 
    % the function measured momentum of the subsequent updates and if the
    % correlation between then is high, it will use this information to
    % accelerate the update in the direction of average velocity 
    
    shift = squeeze(shifts_memory(end,:,:)); 
    
    if ~any((acc_axes & any(shift~=0))); return ; end
    Nmem = size(shifts_memory,1)-1;
    for jj =  find(acc_axes)  % apply only for horizontal, vertical seems to be too unstable 
        if all(shift(:,jj)==0); continue; end
        for ii = (1:Nmem)
            C(ii) = corr(shift(:,jj), squeeze(shifts_memory(ii,:,jj))');
        end
        
        % estimate optimal friction from previous steps 
        decay = fminsearch(  @(x)norm(C - exp(-x*[Nmem:-1:1])), 0); 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        alpha = 2;                          % scaling of the friction , larger == less memory 
        gain = 0.5;                   % smaller -> lower relative speed (less momentum)
        friction = min(1,max(0,alpha*decay));   % smaller -> longer memory, more momentum 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % update velocity map 
        velocity_map(:,jj) = (1-friction)*velocity_map(:,jj) + shift(:,jj);
        % update shift estimates 
        shift(:,jj) = (1-gain) * shift(:,jj)  + gain*velocity_map(:,jj); 
    end

    acc = math.norm2(shift(:,acc_axes)) ./ math.norm2(squeeze(shifts_memory(end,:,acc_axes))); 
    utils.verbose(0,'Momentum acceleration %4.2fx    friction %4.2f', acc, friction )

end


function sinogram = unwrap_data(sinogram, method, boundary)
    % auxiliary function to perform data unwrapping, see
    % math.unwrap2D_fft for detailed help 
    switch lower(method)
        case 'fft_1d'
            % unwrap the data by fft along slices 
            sinogram = -math.unwrap2D_fft(sinogram, 2, boundary);
        case {'none', 'diff'}
            % do nothing 
        otherwise
            error('Missing method')
    end
end

function [shift, err] = find_optimal_shift(sinogram_model, sinogram,weights, MASS, par)
        % given the sinogram_model, measured sinogram, and importance weight for each pixel it tries to
        % estimate the most optimal shift betweem sinogram_model and
        % sinogram to minimize weighted difference || W * (sinogram_model - sinogram + alpha * d(sino)/dX )^2 ||
        
        import math.* 
        shift_x = zeros(size(sinogram_model,3),1,'single');
        shift_y = zeros(size(sinogram_model,3),1,'single');
        
        resid_sino = get_resid_sino(sinogram_model, sinogram, par.high_pass_filter);
        if strcmpi(par.unwrap_data_method, 'none');  resid_sino = imfilter_high_pass_1d(resid_sino,2,par.high_pass_filter,0); end
        smooth_window = 5;
        if par.align_horizontal
            % calculate optimal shift of the 2D projections in horizontal direction 
            dX = get_img_grad_filtered(sinogram_model, 1, par.high_pass_filter, smooth_window );
            if strcmpi(par.unwrap_data_method, 'none');  dX = imfilter_high_pass_1d(dX,2,par.high_pass_filter,0); end
            shift_x = -squeeze(gather(sum2(weights .* dX .* resid_sino) ./ sum2(weights .* dX.^2)));
        end
        clear dX 
        if par.align_vertical
            % calculate optimal shift of the 2D projections in vertical direction 
            dY = get_img_grad_filtered( sinogram_model, 2, par.high_pass_filter,  smooth_window );
            if strcmpi(par.unwrap_data_method, 'none');  dY = imfilter_high_pass_1d(dY,1,par.high_pass_filter,0); end
            shift_y = -squeeze(gather(sum2( weights .* dY .* resid_sino) ./ sum2(weights .* dY.^2)));
        end
        clear dY
        
        shift = [shift_x, shift_y];
       
        if any(isnan(shift(:)))
           warning('Alignment failed, estimated shift is NaN') 
           keyboard
        end
        
        err = squeeze(sqrt(gather(mean2( (weights .* resid_sino).^2)) ))./ MASS; 
end
        

function resid_sino = get_resid_sino(sinogram_model, sinogram, high_pass_filter)
    % calculate filtered difference between sinogram_model  and sinogram 
    % || (sinogram_model  - sinogram) \ast ker ||
    % filtering is used for suppresion of low spatial freq. errors 

    % calculate residuum 
    resid_sino = sinogram_model - sinogram;
    % apply high pass filter => get rid of phase artefacts 
    resid_sino = imfilter_high_pass_1d(resid_sino,2,high_pass_filter); 

end

function d_img = get_img_grad_filtered(img, axis, high_pass_filter,smooth_win)
    % calculate filtered image gradient  between sinogram_model  and sinogram 
    % filtering is used for suppresion of low spatial freq. errors 
    
    import math.* 
    img  = utils.smooth_edges(img,smooth_win, 1+mod(axis,2)); % smooth edges to avoid jumps 
    isReal = isreal(img); 
    Np = size(img);
    if axis == 1
        X = 2i*pi*(fftshift((0:Np(2)-1)/Np(2))-0.5);
        d_img = fft(img,[],2);
        d_img = bsxfun(@times,d_img,X);
        % apply filter in horizontal direction 
        d_img = imfilter_high_pass_1d(d_img,2,high_pass_filter,0,false); 
        d_img = ifft(d_img,[],2);
    end
    if axis == 2
        X = 2i*pi*(fftshift((0:Np(1)-1)/Np(1))-0.5);
        d_img = fft2(img); 
        d_img = bsxfun(@times, d_img, X.');
        % apply filter in horizontal direction 
        d_img = imfilter_high_pass_1d(d_img,2,high_pass_filter,0,false); 
        d_img = ifft2(d_img);
    end
    if isReal; d_img = real(d_img);end 
end

function mask = get_mask(rec, par, circulo)
    % estimate support mask for current reconstruction , assume that there
    % are not holes in the mask 
    
    import utils.* 
    
    [Npix, ~, Nlayers] = size(rec); 
    
    if isempty(par.mask_threshold)
        T = graythresh(rec(:));
    else
        T = par.mask_threshold; 
    end
    
    rec = rec .* circulo; 
    % assume that the sample is roughly vertical pilar -> get only 2D mask
    mask = sum(rec  > T,3) > 0; 
    mask  = imfill(imdilate(mask , strel('disk', 5)),'holes');
    
    % avoid effects of unstrained regions of the reconstruction 
    xt = -Npix/2:Npix/2-1;
    [X,Y] = meshgrid(xt,xt);
    circulo=1-radtap(X,Y,20,round(Npix/2)+20);  
    mask = mask .* circulo; 
    
    plotting.smart_figure(46)
    subplot(2,2,1)
    plotting.imagesc3D(rec , 'init_frame', Nlayers/2)
    colorbar
    axis off image 
    colormap bone 
    title('Current reconstruction')
    subplot(2,2,3)
    plotting.imagesc3D(rec .* ~mask, 'init_frame', Nlayers/2)
    colorbar
    axis off image 
    colormap bone 
    title('Residuum around mask')
    subplot(1,2,2)
    hist(rec(1:100:end), 100)
    plotting.vline(T, '--', 'Current threshold')
    axis tight 
    title('Histogram of the reconstructed values')
    drawnow 

end

function plot_alignment(rec,  sinogram_shifted, weights_shifted, err,shift_upd,shift_total, angles,valid_angles, iter, par)
    % plot results for the current iteration
    % show single sinogram slice, reconstruction slice, evolution of errors
    % and evolution of position correction 
    % ** rec  - reconstructed volume 
    % ** sinogram_shifted - sinogram with applied shifts 
    % ** weights_shifted - importance weights of the sinogram shifted 
    % ** err - error evolution 
    % ** shift_upd - current update of the optimal shift 
    % ** shift_total  - total update of the optimal shifts 
    % ** angles  - angles of the projections 
    % ** valid_angles - bool array of the valid angles, used only to mark the ignored angles in the plot 
    % ** iter - current iteration number 
    % ** par - tomo param structure 

    import utils.*
    import math.*

    [Nlayers,~,Nangles] = size(sinogram_shifted);

    verbose(1,'Plotting')
    fig_id  = 5464; 
    
    if par.windowautopos && ~ishandle(fig_id)  % autopositioning only if the figure does not exists yet 
        plotting.smart_figure(fig_id)
        set(gcf,'units','normalized','outerposition',[0.2 0.2 0.8 0.8])
    else
        plotting.smart_figure(fig_id)
    end
    
    range = gather(math.sp_quantile(rec(:,:,ceil(end/2)), [0.01,0.999], 4)); 
    
    subplot(2,3,1)
    sino_slice = squeeze(sinogram_shifted(ceil(Nlayers/2),:,:))'; 
    if par.is_laminography
        try; sino_slice = sino_slice .* squeeze(weights_shifted(min(ceil(Nlayers/2),end),:,:))'; end
    end
    sino_slice = imfilter_high_pass_1d(sino_slice,2,par.high_pass_filter); 
    
    
    imagesc(sino_slice, math.sp_quantile(sino_slice,[0.01,0.99],5));  
    title(sprintf('High-pass filtered shifted sinogram\nCurrent downsampling: %ix', par.binning'))
    axis off 

    colormap bone 
    
    if par.showsorted
        xaxis = angles; 
        xaxis_label = 'Angle [deg]'; 
    else
        xaxis = 1:Nangles; 
        xaxis_label = '# projection'; 
    end
        

    subplot(2,3,2)
    plot(xaxis, shift_upd(:,1)*par.binning, '.-r')
    hold on 
    plot(xaxis, shift_upd(:,2)*par.binning, '.-b')
    hold off 
    grid on
    legend({'horiz', 'vert'})

    title('Current position update')
    xlim([min(xaxis), max(xaxis)])
    ylabel('Shift  x downsampling [px]')   

    xlabel(xaxis_label)

    subplot(2,3,3)
    plot(xaxis,shift_total(:,1)*par.binning, '.-r')
    hold on 
    plot(xaxis,shift_total(:,2)*par.binning, '.-b')
    hold off 
    title('Total position update')
    legend({'horiz', 'vert'})
    ylabel('Shift x downsampling [px]')
    xlim([min(xaxis), max(xaxis)])
    xlabel(xaxis_label) 
    grid on 

    subplot(2,3,6)
    plot(xaxis(valid_angles), err(iter,valid_angles), 'k.') 
    hold on 
    plot(xaxis(~valid_angles),err(iter,~valid_angles), 'r.') 
    hold off 
    if any(~valid_angles)
        legend({'errors', 'ignored'})
    end
    title('Current error')
    axis tight 
    grid on 
    xlim([min(xaxis), max(xaxis)])
    xlabel(xaxis_label)


    subplot(2,3,5)
    plot(err)
    hold on
    plot(mean(err,2), 'k', 'LineWidth', 3);
    hold off 
    grid on 
    axis tight
    xlim([1,iter+1])
    set(gca, 'xscale', 'log')
    set(gca, 'yscale', 'log')
    title('MSE evolution')
    xlabel('Iteration')
    ylabel('Mean square error')

    subplot(2,3,4)
    Nlayers = size(rec,3);


    if par.is_laminography
        % show also cut in the vertical direction 
        plotting.imagesc3D(rec, 'init_frame', ceil(Nlayers/2), ...
            'fnct', @(x)cat(1, x, ...
            ones(ceil(Nlayers/10), size(rec,2))*range(2), ...
            squeeze(rec(ceil(end/2),:,:))'))
        ylabel('Side view / Top view')
    else
        plotting.imagesc3D(rec, 'init_frame', ceil(Nlayers/2)); 
        ylabel('Horizontal cut')
    end
    axis image
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    caxis(range)
    colormap bone(1024)
    title('Current reconstruction')
    
    drawnow 
   
end
