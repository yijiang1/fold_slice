% INITIALIZE_SOLVER initialize GPU ptycho reconstruction, generate cache values, fftshift data, etc 
% 
% [self, cache] = initialize_solver(self,par) 
% 
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** par       structure containing parameters for the engines 
%
% returns: 
% ++ self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ++ cache     structure with precalculated values to avoid unnecessary overhead

function [self, cache] = init_solver(self,par)
        
    import engines.GPU_MS.shared.*
    import math.*
    import utils.*
    import plotting.*
    import engines.GPU_MS.GPU_wrapper.*
    
    verbose(struct('prefix','GPU/CPU_MS-engine-init'))

    par.Nscans = length(self.reconstruct_ind); %number of scans
    cache.skip_ind = setdiff(1:self.Npos,[self.reconstruct_ind{:}]); %  wrong datasets  to skip 
        
    if ~any(self.probe_support(:))
        self.probe_support = []; 
    end
    %% avoid probe to be larger than a certain oversampling !!!! 
    if isempty(self.probe_support) 
        par.probe_backpropagate = 0;
    end

    if  ~isempty(self.background) && any(self.background(:) > 0)
        Background = self.background;
    elseif  par.background_detection
        Background = 0; 
    else
        Background = [];  %  array of background light 
    end
    
    Noise = [];
    %% prepare data / noise / mask 
    if par.relax_noise &&  ~isempty(self.noise) &&  strcmp(par.likelihood, 'L1')
        Noise = self.noise;
        Noise = (sqrt(posit(self.diffraction + Noise)) - sqrt(posit(self.diffraction - Noise)))/2;
        Noise(self.diffraction == 0) = 1;
        disp('Using measured noise')
        Noise = max(0.5, Noise); 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% PREPARE MASK AND DATA %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    %% prepare mask , note that bool in matlab has size of uint8 !!
    cache.mask_indices = [];
    if any(self.mask(:))
        Mask = [];
        % single mask 
        if all(all(mean(self.mask,3) == self.mask(:,:,1)))
            Mask = self.mask(:,:,1);
        else
            % mask for each scan
            for ll = 1:par.Nscans
                ind = self.reconstruct_ind{ll};
                %if there is only one repeated mask over whole scan 
                if all(all(all(bsxfun(@eq, self.mask(:,:,ind), self.mask(:,:,ind(1))))))
                    Mask(:,:,ll) = self.mask(:,:,ind(1)); 
                end
                cache.mask_indices(ind) = ll;
            end
        end
        if isempty(Mask)
            % mask for each position 
            Mask = self.mask;  % otherwise just store original
            cache.mask_indices(ind) = 1:self.Npos;
        end
        % important to save memory 
        if all(Mask(:) == 1 | Mask(:) == 0)
            Mask = logical(Mask );
        else
            Mask = uint8(Mask*255);  % if there are nonlogical values in mask, store them as uint8 to save memory 
        end
    else
        Mask = [];
    end
    
    %% prepare diffraction data 
    Diffraction = self.diffraction;  % diffraction is intensity, not amplitude, comment by ZC
    if par.upsampling_data_factor
        % downsample the data down to original size to save memory 
        Diffraction = utils.binning_2D(Diffraction, 2^par.upsampling_data_factor) * (2^(2*par.upsampling_data_factor)); 
        if ~isempty(Mask)
            Mask = utils.binning_2D(Mask, 2^par.upsampling_data_factor) == 1; 
        end
    end
    
    Diffraction = single(max(0,Diffraction));

    
    if ~isempty(Mask)
        if size(Mask,3) == par.Nscans && par.Nscans  > 1
            for ll = 1:par.Nscans
                ind = self.reconstruct_ind{ll};
                Diffraction(:,:,ind) = Diffraction(:,:,ind) .* ~Mask(:,:,ll); 
            end
        else
            Diffraction = Diffraction .* ~Mask;
        end
    end
      
    
    if  ~isinf(self.z_distance(end)) %  && mod(Ninf,2)~=0
        % assume inputs already fftshifted, but in case of nearfield
        % fftshift it back for the ASM propagator 
        Noise = fftshift_2D(Noise);
        Diffraction = fftshift_2D(Diffraction);
        Mask = fftshift_2D(Mask);
    end

    
    %%%% compress data if requested %%%%%%
    if par.compress_data
        DATA_MAX = quantile(max2(abs(Diffraction)), 1-1e-2); 
        C_factor_0 = 2;  % compression factor  >=2 seems to be safe, >=4 is pratically lossless 
        if par.compress_data == 1 || DATA_MAX < 2^(2*8) / C_factor_0^2
            Diffraction = sqrt(single(Diffraction));
            if  DATA_MAX < 2^(2*8) / C_factor_0^2
                % simple sqrt compression to 8 bits 
                verbose(1, 'Online data compression to 8-bits')
                Diffraction = uint8(C_factor_0*Diffraction);
                cache.C_factor = C_factor_0;
            elseif DATA_MAX < 2^(2*16) / 16^2
                % failsafe option: sqrt compression to 16 bits 
                verbose(1, 'Online data compression to 16-bits')
                cache.C_factor = 16;  % use compression factor 16, to be super safe just because we have space
                Diffraction = uint16(cache.C_factor*Diffraction);
            else
                error('Online compression will fail')
            end
        elseif par.compress_data == 2
            % SVD subtraction compression to 8 bits (failsafe is compression to 16bits)
            % additionally remove some SVD modes 
            Diffraction = sqrt(single(Diffraction));
            Nmodes = par.Nscans;
            [U,S,V] = fsvd(reshape(Diffraction,prod(self.Np_p),[]), Nmodes);
            ind_relevant = diag(S).^2/sum(diag(S).^2) > 1e-2; % more than 1% of power
            cache.US_diffraction = (U(:,ind_relevant)*S(ind_relevant,ind_relevant));
            cache.V_diffraction = V(:,ind_relevant);
            svd_Diffraction = round(reshape(cache.US_diffraction*cache.V_diffraction',[self.Np_p, self.Npos]));

            %% compress 
            cDiffraction = single(Diffraction) - svd_Diffraction;
            % reestimate optimal compression factor to keep values < 128 
            C_factor = min(C_factor_0, 128/quantile(max2(abs(cDiffraction)), 1-1e-2));
            if C_factor > 3
                verbose(1, 'Online data compression to 8-bits + SVD')
                cache.C_factor = C_factor;
                cache.US_diffraction =  cache.US_diffraction;
                cache.V_diffraction = cache.V_diffraction*C_factor;
                Diffraction = int8(cDiffraction*C_factor);
            elseif DATA_MAX < 2^(2*16) / 16^2
                % sqrt compression to 16 bits 
                %warning(sprintf('Too high online compression of data, it may cause problems\n Compression factor is %2.2f but should be >= 2\n Switching from 8 to 16bits',C_factor))
                verbose(1, 'Online data compression to 16-bits')
                C_factor = 16;
                cache.C_factor = C_factor;
                Diffraction = uint16(C_factor*Diffraction); 
            else
                error('Online compression will fail')
            end
            
            clear svd_Diffraction cDiffraction 

        else
            error('Unimplented level of compression')
        end
    else
       % precalculate sqrt from the data, store as singles 
        Diffraction = sqrt(single(max(0,Diffraction))); % diffraction is amplitude, comment by ZC
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% write back the data arrays 
    self.diffraction = Diffraction; % diffraction is amplitude, comment by ZC
    self.mask = Mask; 
    self.noise = Noise;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% PREPARE GEOMETRY, PROPAGATION, MODES%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % precalculate ASM factor for propagation distance recovery 
    [ASM_difference] = near_field_evolution_gradient(self.Np_p, self.lambda,  self.pixel_size .*self.Np_p );
    cache.ASM_difference = fftshift(ASM_difference);
    
    % custom propagator to account for tilted plane sample
    if any(par.p.sample_rotation_angles(1:2)) && check_option(par.p, 'apply_tilted_plane_correction', 'propagation') 
        % get propagators to the tilted plane 
        [tilted_plane_propagate_fwd, tilted_plane_propagate_back] = ...
            get_tilted_plane_propagators(Garray(self.probe{1}), ...
                                        [par.p.sample_rotation_angles(1:2),0],...
                                        self.lambda, self.pixel_size); 
    else
        tilted_plane_propagate_fwd = [];  tilted_plane_propagate_back = [];
    end
    
  
    if ~iscell(self.affine_matrix)
        self.affine_matrix = {self.affine_matrix}; 
    end
    
    %Note: par.Nmodes is # of probe modes (p.probe_modes). Assigned in
    %load_from_p.m
%     modes = cell(max(par.Nmodes, par.Nlayers),1);
    modes = cell(par.Nlayers,1); % corrected by Zhen Chen
    
    %Comment by YJ: par.Nmodes seems to be # of modes for A-fly scan
%     for i = 1:max(par.Nmodes, par.Nlayers) 
    %par.Nlayers equals to the # of slices in the object (excluding the
    %vacuum layer)
    for i = 1:par.Nlayers % modified by ZC
        verbose(2,'Creating new modes files ')
        modes{i}.lambda =  self.lambda;
        
        % decompose affine matrix into scale, asymmetry, rotation, shear 
        %affine = scale*[1+asym/2,0; 0,1-asym/2]*[cosd(rot), sind(rot); -sind(rot), cosd(rot)] * [1,0;tand(shear),1];
        
        affine_matrix = self.affine_matrix{min(i,end)};
        [scale, asymmetry, rotation, shear]  = decompose_affine_matrix(affine_matrix);
        
        % store initial geometry parameters 
        modes{i}.scales =  repmat(scale, 1,par.Nscans); 
        modes{i}.asymmetry =  repmat(asymmetry, 1,par.Nscans); 
        modes{i}.shear =   repmat(shear, 1,par.Nscans); 
        modes{i}.rotation =  repmat(rotation, 1,par.Nscans); 
        modes{i}.affine_matrix =  repmat(affine_matrix, 1,1,par.Nscans); 
        modes{i}.shift_scans =  zeros(2, par.Nscans); 
        modes{i}.probe_scale_upd = 0; 
        modes{i}.probe_rotation = ones(1,par.Nscans) * par.sample_rotation_angles(3); % one rotation per scan 
        if par.mirror_objects
            modes{i}.probe_rotation = modes{i}.probe_rotation .* [1,-1];  % flip the coordinates for mirrored object (ie 0 vs 180deg rotation)
        end
        modes{i}.probe_rotation_all = zeros(self.Npos,1);
        for jj = 1:par.Nscans
            modes{i}.probe_rotation_all(self.reconstruct_ind{jj}) = modes{i}.probe_rotation(jj); % one rotation per scan 
        end
        
        distance = self.z_distance(min(end,i));
        
        if ~isinf(distance)
            verbose(2, 'Layer %i distance %g um  ', i, distance*1e6 )
        end
        modes{i}.distances = distance;
        if is_used(par, 'fly_scan') && (~isfield(modes{i}, 'probe_positions') || isempty(modes{i}.probe_positions) )
            %% get positions for fly scans 
             self = prepare_flyscan_positions(self, par); 
             modes{i}.probe_positions = self.modes{i}.probe_positions; %added by YJ. seems like a bug
             modes{i}.probe_positions_0 = self.probe_positions_0; %added by YJ. seems like a bug
        else
           %% get positions for normal tomo 
            try  %  try to reuse the positions of there are saved 
                modes{i}.probe_positions = self.modes{i}.probe_positions;
                verbose(2,'Using saved positions')
            catch
                if (modes{i}.scales(end) == modes{1}.scales(end)) && ~isempty(self.probe_positions)
                    modes{i}.probe_positions = self.probe_positions;
                    verbose(0,'Using saved positions')
                else
                    verbose(2,'Using original positions')
                    modes{i}.probe_positions = (affine_matrix*self.probe_positions_0')';
                end
            end
            try
                modes{i}.probe_positions_0 = self.modes{i}.probe_positions_0;
            catch
                modes{i}.probe_positions_0 = self.probe_positions_0;
            end

        end
        modes{i}.probe_positions_update = { zeros(size(modes{i}.probe_positions)) };
        modes{i}.probe_positions_all = {modes{i}.probe_positions};
        modes{i}.probe_positions_weight = zeros(self.Npos, 1);   
        if isfield(self, 'probe_fourier_shift') && ~isempty(self.probe_fourier_shift)  && i == 1
            modes{i}.probe_fourier_shift = self.probe_fourier_shift; 
        else
            modes{i}.probe_fourier_shift =  zeros(self.Npos,2);
        end

        if ~isempty(self.probe_support) && i <= par.Nrec
            modes{i}.probe_support = self.probe_support; 
            if i == 1
                verbose(2,'Using real-space probe support')
            end
        else 
            modes{i}.probe_support = [];
        end
        
        if ~isempty(self.probe_support_fft) && i <= par.Nrec && ~check_option(par,'probe_support_tem')
            modes{i}.probe_support_fft = fftshift(self.probe_support_fft); 
            if i == 1
                verbose(2,'Using far-field probe support')
            end                      
        elseif check_option(par,'probe_support_tem') % not shift for TEM aperture, by Zhen Chen
            modes{i}.probe_support_fft = self.probe_support_fft; 
        else 
            modes{i}.probe_support_fft = [];
        end
        
        F = mean( self.pixel_size)^2 .* mean(self.Np_p) /   (modes{i}.lambda * modes{i}.distances); 
        if F ~= 0
            verbose(3,'Nearfield propagation: Fresnel number/Npix %3.3g', F)
        end
        scale = modes{i}.scales(end);
        modes{i}.ASM_factor = [] ;
        modes{i}.cASM_factor = [] ;

        if ~isinf(modes{i}.distances(end)) % Forward Fresnel propagator in k-space, ASM, commented by ZC
            %% near field factor            
            %ASM =  exp( modes{i}.distances(end)* cache.ASM_difference);
            % modified by YJ: use H instead of dH (which is an approximation)
            [~,ASM,~,~] = near_field_evolution(ones(self.Np_p), modes{i}.distances(end), self.lambda,  self.pixel_size .*self.Np_p, true );
            ASM = fftshift(ASM);
            modes{i}.ASM_factor = ASM;
            modes{i}.cASM_factor = conj(ASM);
        end
        
        %% far field factor
        modes{i}.FAR_factor = [];
        modes{i}.cFAR_factor = conj(modes{i}.FAR_factor);
  
        if isinf( par.probe_backpropagate)
             modes{i}.support_fwd_propagation_factor = inf; 
             modes{i}.support_back_propagation_factor = -inf;
        elseif par.probe_backpropagate ~= 0
            [~, modes{i}.support_fwd_propagation_factor] = utils.prop_free_nf( self.probe{1}(:,:,1), par.probe_backpropagate,...
                modes{i}.lambda,  self.pixel_size ./ scale );
  
            modes{i}.support_fwd_propagation_factor = fftshift( modes{i}.support_fwd_propagation_factor );
            modes{i}.support_back_propagation_factor = conj(modes{i}.support_fwd_propagation_factor);
        else
            modes{i}.support_fwd_propagation_factor = [];
            modes{i}.support_back_propagation_factor = [];
        end
        
        modes{i}.tilted_plane_propagate_fwd = tilted_plane_propagate_fwd;
        modes{i}.tilted_plane_propagate_back = tilted_plane_propagate_back;
    end
    
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% PREPARE PROBES, INCOHERENT MODES %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    probe_0 = mean(self.probe{1},3);
    probe = cell(1,par.probe_modes);  % Added by ZC and fixed by YJ
    for i = 1:par.probe_modes
        try
            probe{i} = self.probe{i};
                %% test if the probe size ok for the variable probe settings etc 
            assert( size(probe{i},4) == 1+par.variable_probe_modes || ...
                ~(par.variable_probe) || i > 1)
            assert((size(probe{i},3) ==1 || par.variable_probe) || ...
                   (size(probe{i},3) == par.Nscans && ~par.share_probe)  )  % no variable prob extension and multiple probes used 
            assert(size(probe{i},3) == par.Nscans || par.share_probe || par.variable_probe, 'Wrong probe size for not shared probe option')
        catch
            if i <= par.Nrec || is_used(par, 'fly_scan')
                verbose(2, 'Creating probe')

                if ~par.share_probe && size(probe{i},3) == 1
                     % dont share probe between scans 
                     probe{i} = repmat(probe_0,[1,1,par.Nscans]);
                end
                if  (par.variable_probe && par.variable_probe_modes > 0) && i == 1
                    verbose(2,'Creating variable probe ')
                    probe{i}(:,:,:,2:1+par.variable_probe_modes) = ...
                        randn([self.Np_p, size(probe{i},3), par.variable_probe_modes])+randn([self.Np_p,size(probe{i},3), par.variable_probe_modes])*1i;
                    continue 
                end
            end
            if length(probe) < i   % none of above 
               % simply create slightly shifted modes in fourier domain, it is useful for
               % inital guess of incoherent modes after orthogonalization 
                step = median(diff(self.probe_positions_0)); 
                probe{i} = 0.01*fftshift(imshift_fft(fftshift(probe_0), randn, randn, false));
            end
            
            % fill the unreconstructed positions if the OPRP method is used
            if par.variable_probe && is_method(par, 'PIE') && i ==1 
                ind_wrong = setdiff(1:self.Npos, [self.reconstruct_ind{:}]);
                probe{i}(:,:,ind_wrong) = repmat(mean(probe{i},3),1,1,length(ind_wrong));
            end
         end
    end
    if par.probe_modes > par.Nrec
        %  orthogonalization of incoherent probe modes 
        if is_used(par, 'fly_scan')
            probe_tmp = probe;
            % orthogonalize the modes with all the other shifted modes 
            for i = 1:par.Nrec
                dx = median(modes{i}.probe_positions - modes{1}.probe_positions);
                probe_tmp{i} = imshift_fft(probe_tmp{i}, dx); 
            end
            probe_tmp = ortho_modes(probe_tmp);  % perform othogonalization 
            probe(1+par.Nrec:par.probe_modes) = probe_tmp(1+par.Nrec:par.probe_modes);
        else
            ind = [1,1+par.Nrec:par.probe_modes];  % skip polyvave/multilayer probe_tmp 
            probe(ind) = ortho_modes_eig(probe(ind));  %% slightly better  
        end
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% PREPARE OBJECT, MULTILAYER OBJECT %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    %% updated illumination
    aprobe2 = abs(self.probe{1}(:,:,1)).^2; 
    for ll = 1:par.Nscans
        if par.share_object
            ind = [self.reconstruct_ind{:}];
        else
            ind = self.reconstruct_ind{ll};
        end
        [cache.oROI_s{1}] = find_reconstruction_ROI( modes{1}.probe_positions,self.Np_o, self.Np_p);
        % avoid oscilations by adding momentum term 
        illum_sum_0{ll} = Ggather(set_views(Gzeros(self.Np_o), Garray(aprobe2), 1,1, ind, cache));
    end
    
    %% multilayer extension 
    % modified by YJ for more dynamic initialization

    % Step 1: choose specific layers from initial object file
    N_layer_input_obj = size(self.object,2);
    par.init_layer_select(par.init_layer_select<0) = [];
    par.init_layer_select(par.init_layer_select>N_layer_input_obj) = [];
    if ~isempty(par.init_layer_select)
        object_temp = cell(size(self.object,1),length(par.init_layer_select));
        for ll = 1:par.Nscans
            for jj=1:length(par.init_layer_select)
                object_temp{ll,jj} = self.object{ll,par.init_layer_select(jj)};
            end
        end
        self.object = object_temp;
    end
    
    % Step 2: pre-process layers
    switch par.init_layer_preprocess
        case 'avg' % only use the averaged layer
            verbose(0,'Average initial layers')
            for ll = 1:par.Nscans
                obj_avg = prod(cat(3,self.object{ll,:}),3); 
                obj_avg = abs(obj_avg).*exp(1i*phase_unwrap(angle(obj_avg))/size(self.object,2));
                for jj=1:size(self.object,2)
                    self.object{ll,jj} = obj_avg;
                end
            end
        case 'interp' % interpolate layers 
            if ~isempty(par.init_layer_interp)
                verbose(0,'Interpolate %d initial layers to %d layers', size(self.object,2), length(par.init_layer_interp))
                for ll = 1:par.Nscans
                    obj_temp = cat(3,self.object{ll,:});
                    [N_obj_y,N_obj_x,N_obj_z] = size(obj_temp); 
                    [X,Y,Z] = meshgrid(linspace(1,N_obj_x,N_obj_x),linspace(1,N_obj_y,N_obj_y),linspace(1,N_obj_z,N_obj_z));
                    [Xq,Yq,Zq] = meshgrid(linspace(1,N_obj_x,N_obj_x),linspace(1,N_obj_y,N_obj_y),par.init_layer_interp);
                    obj_temp = interp3(X,Y,Z,obj_temp,Xq,Yq,Zq,'spline');
                    for jj=1:size(obj_temp,3)
                        self.object{ll,jj} = obj_temp(:,:,jj);
                    end
                end
            end
        case {'','all'} % default: keep all layers
            % nothing to do
        otherwise
            error('Invalid init_layer_preprocess!')
    end
    
    % Step 3: add or remove layers based on par.Nlayers
    if size(self.object,2) > par.Nlayers
        warning('Initial object has more layers than Nlayers')
        for ll = 1:par.Nscans
            self.object{ll,1} = prod(cat(3,self.object{ll,:}),3); 
        end
        self.object(:,2:end) = [];
    end
    
    if size(self.object,2) < par.Nlayers
        N_add = par.Nlayers - size(self.object,2);
        verbose(0,'Add %d more layers from %d layer(s)', N_add, size(self.object,2))
        for ll = 1:size(self.object,1) %loop over scans
            obj{ll} = self.object(ll,:);
            switch par.init_layer_append_mode
                case 'avg' %Not sure when this is useful, but I'll keep it for now
                    verbose(0,'Append averaged layer')
                    obj_pre = prod(cat(3,self.object{ll,:}),3); 
                    obj_pre = abs(obj_pre).*exp(1i*phase_unwrap(angle(obj_pre))/size(self.object,2));
                    obj_post = obj_pre;
                case 'edge'
                    verbose(0,'Append 1st/last layer')
                    obj_pre = self.object{ll,1};
                    obj_post = self.object{ll,end};
                case {'','vac'}
                    verbose(0,'Append vacuum layer')
                    %obj_pre = ones(self.Np_o, 'single') + 1e-9i*randn(self.Np_o, 'single');
                    obj_pre = ones(self.Np_o, 'single');
                    obj_post = obj_pre;
                otherwise
                    error('Invalid init_layer_append_mode!')
            end
            for ii = 1:N_add
                if mod(ii, 2) == 1
                    obj{ll}{end+1} = obj_post; % add slice at the end 
                else
                    obj{ll}(2:end+1) = obj{ll};
                    obj{ll}{1} = obj_pre; % add slice at the beginning 
                end
            end
        end
        self.object = cat(1, obj{:}); %combine all scans
    end
    
    % if object has more layers but only one is needed
    if size(self.object,2) > 1 && par.Nlayers == 1
       for ll = 1:par.Nscans
           object{ll,1} = prod(cat(3,self.object{ll,:}),3); 
       end
       self.object = object; 
    end

    % At this point: size(self.object,3) should equal to par.Nlayers
    %{ 
    %% MO's code, should be useless now. I'll keep it for now in case of
    bugs in steps 1-3. 
    for j = 1:par.Nlayers  
        for i = 1:max(1, par.Nscans * ~par.share_object) % loop over scans
            try
                object{i,j} = self.object{min(end,i),j};
                object{i,j}(1);
            catch
                verbose(0, 'add transparent slice') % add extra layers 
                %% add fully transparent slice at the end 
                object{i,j} = ones(self.Np_o, 'single');
                if size(self.object,2) == 1
                    % swap order of the new layers to keep the original
                    % reconstruction in center
                    object(i,:) = object(i,end:-1:1); 
                end
            end
        end
    end
    %}
    % Step 4: assign self.object to object and rescale layers if needed
    if par.init_layer_scaling_factor~=1
        verbose(0,'Rescale each layer by %f', par.init_layer_scaling_factor)
    end
    for j = 1:par.Nlayers
        for i = 1:max(1, par.Nscans * ~par.share_object) % loop over scans
            if par.init_layer_scaling_factor~=1
                object_temp = self.object{min(end,i),j};
                object_temp_ph = phase_unwrap(angle(object_temp))*par.init_layer_scaling_factor;
                object{i,j} = abs(object_temp).*exp(1i.*object_temp_ph);
            else
                object{i,j} = self.object{min(end,i),j};
            end
        end
    end
    
    for i = 1:numel(object)
        object{i} = single(object{i});
        object{i} =  complex(object{i});
    end

    for i = 1:numel(probe)
        probe{i} = single(probe{i});
        probe{i}  = complex(probe{i});
    end
    
    %% STORE RESULTS TO SELF CLASS
    self.object = object; 
    self.probe = probe; 
    self.modes = modes; 
    self.diffraction = Diffraction; 
    self.noise = Noise; 
    self.mask = Mask; 
    self.background = reshape(Background,1,1,[]);        

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% PRECALCULATE USEFUL VALUES %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
    if ~isfield(self, 'probe_evolution' ) 
        % initial coefficients for OPRP approximation 
        self.probe_evolution(:,1) = ones(self.Npos,1);  % first mode is constant 
    end
    new_probe_modes_ind = 1+(size(self.probe_evolution,2):par.variable_probe_modes); 
    self.probe_evolution(:,new_probe_modes_ind) = 1e-6*randn(self.Npos,length(new_probe_modes_ind));
           
    if par.variable_probe 
        pnorm = norm2(self.probe{1});
        self.probe{1}(:,:,:,2:end) = self.probe{1}(:,:,:,2:end) ./ pnorm(1,1,:,2:end);
        self.probe_evolution(:,2:end) = self.probe_evolution(:,2:end) .* squeeze(mean(pnorm(1,1,:,2:end),3))';
    end
    

    if par.background_detection || ~isempty(self.background)
        %% auto-estimate background correction distribution 
        if isempty(self.mask)
            mask = 0;
        else
            mask = self.mask; 
        end
        if par.background_detection
            background_weight = sum(self.diffraction.^2,3) ./ max(1,sum(~mask,3)); 
            background_weight = imgaussfilt(background_weight,1);
            background_weight = 1./sqrt(max(1e-3, background_weight)) .* ~any(mask,3);
            background_weight = max(0,background_weight - 0.3*mean(background_weight(:)));
            cache.background_weight = ( fftshift_2D(background_weight / sum2(background_weight)));
        end
        
        if isinf(par.background_width)
            cache.background_profile = 1;
        else 
            mdiffr = Garray(fftshift(mean(get_modulus(self,cache,1:self.Npos,false).^2,3)));

            W = par.background_width;
            X = (-self.Np_p(1):self.Np_p(1)-1);
            Y = (-self.Np_p(2):self.Np_p(2)-1);
            [X,Y] = meshgrid(X,Y);

            background_profile = exp(-sqrt( (X/W(1)).^2 +(Y/W(1)).^2)); 
            background_profile = conv2(mdiffr,background_profile, 'same');
            background_profile = background_profile / max2(background_profile);
            
            background_profile = utils.crop_pad(background_profile,self.Np_p); 
            cache.background_profile = gather(fftshift(background_profile)); 
 
        end
         
        if ~isempty(self.diffraction_deform_matrix)
            apply_deform = @(x,D)single(reshape(full(D * double(reshape(x,[],size(x,3)))), size(x)));
            % apply deformation effects caused by tilted sample , be sure to enforce the mask before
            % interpolation, the hotpixels can spread around after the correction
            if isscalar(cache.background_profile)
                cache.background_profile = ones(self.Np_p, 'single'); 
            end
            cache.background_profile = apply_deform(cache.background_profile, self.diffraction_deform_matrix');
        end
    else
         cache.background_profile_weight = 1; 
    end

    for  ll = 1:par.Nscans
        illum_sum_0{ll} = Ggather(illum_sum_0{ll}); 
        cache.MAX_ILLUM(ll) = max(illum_sum_0{ll}(:));
        cache.illum_sum_0{ll} = illum_sum_0{ll}; 
    end

    %% precalculate illumination ROIs
    cache = precalculate_ROI(self,cache, Ggather(sqrt(aprobe2)));

    %% prepare mask needed for subpixel shifts of object views 
    cache.apodwin = single(0.1+0.9*tukeywin(self.Np_p(1),0.05) .* tukeywin(self.Np_p(2), 0.05)');

    if par.initial_probe_rescaling
        %% initial rescaling of probe intensity , just a very rough guess     
        % modified by ZC, propagate probe to far field
        mean_aPsi = mean2(abs(fft2_safe(self.probe{1}(:,:,1))).^2); 
        % old method: self.modes{end} is the vaccum layer
        %mean_aPsi = mean2(abs(fwd_fourier_proj(self.probe{1}(:,:,1), self.modes{end})).^2); 
        
        mean_diffraction_intensity = mean(mean2(self.diffraction(:,:,randi(self.Npos, [10,1])).^2));  % take roughly average intensity  % bug fixed by Zhen Chen, previous no ^2

        for ii = 1:par.probe_modes
            self.probe{ii} = self.probe{ii} * sqrt( mean_diffraction_intensity / mean_aPsi); 
        end   
    end

end
