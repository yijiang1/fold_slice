function [ p, fdb ] = ML_MS_rot( p )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

    global opt_time
    fdb.status = core.engine_status;
    
    opt_time = 0;
    verbose(1, 'Starting multi-slice non-linear optimization')
    
    % ===== 3ML =====
    N_layer = p.N_layer;
    Ny = p.asize(1);
    Nx = p.asize(2);
    lambda = p.lambda; 
    k  = 2*pi/lambda; 
    
    % ----- Calculate the propagator
    [Xp,Yp] = meshgrid(([1:p.asize(2)]-floor(p.asize(2)/2)+1)*p.dx_spec(2), ([1:p.asize(1)]-floor(p.asize(1)/2)+1)*p.dx_spec(1));
    Xp = ifftshift(Xp); 
    Yp = ifftshift(Yp);
    dx = Xp(1,2)-Xp(1,1);
    Fx = Xp/(Nx*dx^2);
    dy = Yp(2,1)-Yp(1,1);
    Fy = Yp/(Ny*dy^2);
%     for n = 1:N_layer-1
%         propagation{n} = exp( 1j*k*p.delta_z(n)*sqrt( 1-(lambda*Fx).^2-(lambda*Fy).^2 ) );
%         propagation_back{n} = exp( 1j*k*(-p.delta_z(n))*sqrt( 1-(lambda*Fx).^2-(lambda*Fy).^2 ) );
%     end
    p.Fx = Fx;
    p.Fy = Fy;
    
    core.errorplot; % clear the persistent variable
    for outer = 1:p.ms_outer_iter
    fprintf('##### Outer %d #####\n', outer);
        
    % ----- Initialization
    if (length(p.ms_init_ob_fraction) ~= p. N_layer) || sum(p.ms_init_ob_fraction)~=1
        fprintf('-- (Initialization) p.ms_init_ob_fraction bad, will use 1/N_layer for all layers \n');
        p.ms_init_ob_fraction = ones(1,p. N_layer)/p. N_layer;
    end
           
    for obnum = 1:p.numobjs
        if strcmp(p.initial_iterate,'file') || var(angle(p.object{obnum}(:)))>0.1 % spicify input or propagating results from another engine
            ob_phase{obnum} = engines.ML_MS.fun_ramp_unwrap(p.object{obnum}, p.asize);
        else
            ob_phase{obnum} = angle(p.object{obnum});
        end
        for n = 1:N_layer            
            for obmode = 1:p.object_modes            
                object_layer{obnum}{obmode}{n} = abs(p.object{obnum}).^(p.ms_init_ob_fraction(n)) .* exp(1i*ob_phase{obnum}.*p.ms_init_ob_fraction(n));
            end
        
            if p.use_display
                figure(99+obnum); subplot(N_layer,2,2*n-1); imagesc(abs(object_layer{obnum}{1}{n})); colormap bone; axis equal xy tight; title(['amplitude, layer' num2str(n)]); caxis([0 2]); colorbar; drawnow;
                figure(99+obnum); subplot(N_layer,2,2*n); imagesc(angle(object_layer{obnum}{1}{n})); colormap bone; axis equal xy tight; title('phase'); caxis([-pi pi]); colorbar; drawnow;
            end
        end
    end	
    probes = p.probes;
        
    recon_time_tic = tic;
    recon_time = [];
    delta_z_iter = [];   
    % ========== 
    
    if p.probe_mask_bool
        
        if p.probe_mask_use_auto
            verbose(2, 'Using a probe mask from probe autocorrelation.');
            to_threshold = -real(auto);
        else
            verbose(2, 'Using a circular probe mask.');
            [x,y] = meshgrid(-p.asize(2)/2:floor((p.asize(2)-1)/2),-p.asize(1)/2:floor((p.asize(1)-1)/2));
            to_threshold = (x.^2 + y.^2);
            clear x y
        end
        to_threshold_flat = reshape(to_threshold, [prod(p.asize) 1]);
        [~, ind] = sort(to_threshold_flat);
        probe_mask_flat = zeros([prod(p.asize) 1]);
        probe_mask_flat(ind(1:ceil(p.probe_mask_area * prod(p.asize)))) = 1;
        p.probe_mask = reshape(probe_mask_flat, p.asize);
        clear to_threshold to_threshold_flat dummy ind probe_mask_flat
    else
        p.probe_mask = ones(p.asize);
    end
    
    % Taking care to pass some needed functions in p
    fnorm = sqrt(prod(p.asize));

    %%% Optimization error metric
    if isfield(p,'opt_errmetric'),
        switch lower(p.opt_errmetric)
            case 'l1'
                verbose(1, 'Using ML-L1 error metric'),
            case 'l2'
                verbose(1,'Using ML-L2 error metric'),
            case 'poisson'
                verbose(1,'Using ML-Poisson'),
            otherwise
                error([p.opt_errmetric ' is not defined'])
                return;
        end
    else
        p.opt_errmetric = 'poisson';
        verbose(1, 'Using default Poisson error metric')
    end
    
    %%% Set specific variables needed for different metrics %%%
    switch lower(p.opt_errmetric)
        case 'poisson'
            fmag2 = p.fmag.^2;
            fmag2renorm = fmag2/p.renorm^2;
            initialerror = p.renorm^2*sum( p.fmask(:).*( (fmag2renorm(:)+0.5).*log(fmag2renorm(:)+1) ...
                - fmag2renorm(:) ...
                - 1 + 0.5*log(2*pi) + 1./(12*(fmag2renorm(:)+1))  ...
                - 1./(360*(fmag2renorm(:)+1).^3) + 1./(1260*(fmag2renorm(:)+1).^5) )) ... %% Approximation to log(n!) http://www.johndcook.com/blog/2010/08/16/how-to-compute-log-factorial/
                + sum( fmag2(:)*log(renorm^2) );
            clear fmag2renorm
        case 'l1'
            initialerror = 0;
            fmag2 = 0;
        case 'l2'
            initialerror = 0;
            fmag2 = p.fmag.^2;
        otherwise
            error(['Error metric ' p.opt_errmetric 'is not defined'])
    end
    
    
    %%% Regularization
    Npix = 0;
    if p. reg_mu > 0
        for obnum = 1:p.numobjs
            Npix = Npix + p.object_size(obnum,1)*p.object_size(obnum,2);
        end
        Nm = prod(p.asize)*size(p.fmag,3);
        K = 8*Npix^2/(Nm*p.Nphot);
        creg = p.renorm^2*p.reg_mu/K;
    else
        creg = 0;
    end
    
    %%% Sieves preconditioning
    if any(p.smooth_gradient) ~= 0
        if length(p.smooth_gradient) <= 1
            % Hanning regularization
            auxi = fract_hanning_pad(512,512,0);
            auxi = fftshift(ifft2(auxi));
            smooth_gradient = real(auxi(256:258,256:258));  % Regularization kernel ( = 0 to omit)
        end
    else
        smooth_gradient = 0;
    end
    
    %% ===== 3ML main =====
    p.ms_opt_flags_local = p.ms_opt_flags;  
    
    if p.ms_opt_flags(3)
        N_iter_outer = ceil(p.ms_opt_iter/p.ms_opt_z_param(1)) + floor(p.ms_opt_iter/p.ms_opt_z_param(1));
    else
        N_iter_outer = 1;
        N_iter_inner = p.ms_opt_iter;
    end
    
    %core.errorplot; % clear the persistent variable
    for iter_outer = 1:N_iter_outer
        
        if p.ms_opt_flags(3)         
            p.ms_opt_flags_local(3) = ~mod(iter_outer,2); % Alternates between updating delta_z, [0 1 0 1...]
            if p.ms_opt_flags_local(3)
                N_iter_inner = p. ms_opt_z_param(2); % when update delta_z
            else
                N_iter_inner = p. ms_opt_z_param(1); 
            end
        end
       
        optimize_object_layer = p.ms_opt_flags_local(1);
        optimize_probes = p.ms_opt_flags_local(2);
        optimize_delta_z = p.ms_opt_flags_local(3);
        delta_z_iter = [delta_z_iter; p.delta_z(:)'];
        
        % -- Arranging optimization vector 
        xopt = [];
        if optimize_object_layer
            for obnum = 1:p.numobjs
                for obmode = 1:p.object_modes
                    for n = 1:N_layer
                        xopt = [xopt; reshape([real(object_layer{obnum}{obmode}{n}(:)).'; imag(object_layer{obnum}{obmode}{n}(:)).'], [], 1)];
                    end
                end
            end
        else
            p.object_layer = object_layer;
        end
        if optimize_probes
            xopt = [xopt; reshape([real(probes(:)).'; imag(probes(:)).'], [], 1)];
        else
            p.probes = probes;
        end
        if optimize_delta_z
            xopt = [xopt; p.delta_z];
        end   
      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Main optimization loop %%%
        opt_time = tic;
        [tmp, p] = engines.ML.cgmin1('engines.ML_MS.gradient_ptycho_MS', xopt, N_iter_inner, p.opt_ftol, p.opt_xtol,...
            p, p.fmag, fmag2, p.fmask, p.numobjs, p.object_size, p.numprobs,...
            p.numscans, p.scanindexrange, initialerror, fnorm, p.probe_mask, p.plot.errtitlestring, p.plot_mask,...
            p.plot_ind, creg, smooth_gradient);

        opt_time = toc(opt_time);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
   
        % --- Error and time
        ms_opt_error = core.errorplot([]);   % Only reads the persistent variable
        ms_error(2, iter_outer) = ms_opt_error(end);
        recon_time(iter_outer) = toc(recon_time_tic);
        
        fprintf('[iter_outer %d, with flags %d %d %d] opt_time = %.2f min, total %.0f min\n\n', ...
            iter_outer, p.ms_opt_flags_local(1), p.ms_opt_flags_local(2), p.ms_opt_flags_local(3), opt_time/60, recon_time(end)/60);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Arrange solution vector %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if optimize_object_layer
            for obnum = 1:p.numobjs
                for obmode = 1:p.object_modes
                    for n = 1:N_layer
                        o_size = p.object_size(obnum, :);
                        o_numel = prod(o_size);
                        object_layer{obnum}{obmode}{n} = reshape(tmp(1:2:2*o_numel), o_size) + ...
                            1i*reshape(tmp(2:2:2*o_numel), o_size);
                        tmp = tmp(2*o_numel+1:end);
                    end
                end
            end
        end
        if optimize_probes
            probeelements = [p.asize p.numprobs p.probe_modes];
            probes = reshape(tmp(1:2:2*prod(probeelements)), probeelements) + ...
                1i*reshape(tmp(2:2:2*prod(probeelements)), probeelements);
            tmp = tmp(2*prod(probeelements)+1:end);
        end    
        if optimize_delta_z
            delta_z = tmp;
            tmp = tmp(size(delta_z)+1:end);
            fprintf(' Optimized delta_z = %.4f um\n',delta_z(1)*1e6);
            if delta_z<0
                fprintf(' Force delta_z = 0 um\n');
                p.delta_z = 0;
            else
                p.delta_z = delta_z;
            end
        end  
        
        if ~isempty(tmp)
            warning('Temporary vector is not empty, optimized values not assigned');
        end
    
    end % end iter_outer
        
    verbose(2, 'Finished');
    verbose(2, 'Time elapsed in optimization refinement: %f seconds', opt_time);
    
    p.delta_z_iter = delta_z_iter;
    p.error_metric.iteration = 1:length(ms_opt_error); 
    p.error_metric.value = ms_opt_error; 
    p.error_metric.err_metric = p.opt_errmetric; 
    p.error_metric.method = p.name;
    p.recon_time = recon_time;
    
    delta = (p.meta{2}.spec.fsamroy - p.meta{1}.spec.fsamroy)/2;
    for ii = 1:p.numscans
        if p.share_object
            obnum = 1;
        else
            obnum = ii;
        end
        if p.share_probe
            prnum = 1;
        else
            prnum = ii;
        end
        object = ones(size(object_layer{obnum}{:}{n})); 
        for n = 1:N_layer 
            object = object .* object_layer{obnum}{:}{n};   % Combine all layers
            p.object_layers{obnum}{n} = object_layer{obnum}{:}{n}; % For each scan (object number) 
        end   
        p.object{obnum} = object;
        p.probes = probes(:,:,prnum,:);
    end
    for obnum = 1:p.numobjs            
        p.proj_rot{obnum} = MS.fun_generate_prec_proj_shift(p, obnum, delta*(-1)^(obnum-1), [1 0]); % rotate to the middle angle
    end
    [p.proj_rot{1}, p.proj_rot{2}, delta_all] = fun_align_img(p.proj_rot{1}, p.proj_rot{2}, p.asize, p.dx_spec(1)); % align shifted-objects
    p.object_rot =  MS.fun_average(p.proj_rot{1}, p.proj_rot{2}, p.asize); % combine to the middle angle
    x = abs(p.proj_rot{2} - p.proj_rot{1}); 
    p.diff_proj(outer) = sum(sum(x));
    fprintf('----- Sum(abs(difference)) = %.2e ----- \n', p.diff_proj(outer))
    
    if outer < p.ms_outer_iter
        temp_obj = MS.fun_generate_prec_proj_shift(p, 2, delta*(-2), [1 0]); % combine from the other angle
        temp_obj = shiftpp2(temp_obj, delta_all(1), delta_all(2));
        p.object{1} = MS.fun_average(p.object{1}, temp_obj, p.asize);

        temp_obj = MS.fun_generate_prec_proj_shift(p, 1, delta*(2), [1 0]); % combine from the other angle
        p.object{2} = MS.fun_average(p.object{2}, temp_obj, p.asize);

        if 0
            img_plot =  p.object{1} - p.object{2};
            figure(21); 
            imagesc(abs(img_plot)); colormap bone; axis xy tight equal; colorbar
        end
    end
             
    end 
    
    %%%%%%%%%%%%%%%%%
    %%% Last plot %%%
    %%%%%%%%%%%%%%%%%
    if p.use_display||p.store_images
        p.plot.extratitlestring = sprintf(' (%dx%d) - ML', p.asize(2), p.asize(1));
        core.analysis.plot_results(p, p.use_display, p.store_images);
    end
    core.errorplot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% end optimization refinement %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



