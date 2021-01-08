%% TEST TEMPLATE FOR FUNTIONALITY CHECK OF MULTILAYER EXTENSION IN GPU ENGINES 
% 1) call standard template to get fresh settings defaults 
% 2) generate artificial data that should serve as a standart test "sample"
% 3) call GPU engine with different basic functionalities and test if all still works


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% IDEA BEHIND THE ITERATIVE PTYCHO TOMOGRAPHY %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Why is it useful ? 
% ==================
% 1) resolution improvement for ideal noise-limited datasets 
% Ptychotomography is gaining information from overlap of the projection in
% the 3D Fourier space of the sample volume. So theoretically, the provided
% constrait s better than for simple 2D overlap in ptychography. 
% But if the collected anglular sampling is following the Crowther criterion,
% there is no overlap between the projection for the highest spatial
% frequencies !!! 
%
% Conclusion: Iterative ptychotomo should not improve resolution for ideal
% noise-limited datasets. It can theoretically improve SNR for the middle
% and low spatial frequencies, but since the error in the low spatial 
% frequencies provided by ptychography is anyway rather limited by
% systematic errors, there is not much to expect. 
% 
% 2) position refinement 
% =======================
% Position errors in the projection, e.g. from drifts, will result in
% errors that are effecting also lower spatial frequencies. Conventional
% ptychography can refine the positions and geometry errors, but it is very
% slow and not too reliable. Much better is to use information from
% neighbouring projections or 180 deg mirrored projections 
% 
% Conclusion: 
% Iterative ptychotomo should help with position refinement and improve
% reconstrution of samples, where this could be limitting, e.g. large catalyst
% particles. 
%
% 3) extended depth of focus 
% ===========================
% Extended depth of focus can be achieve via ptychography but
% reconstruction of N layers is opening N-times more degrees of freedom
% that reconstruction of single layer. This means that eDOF ptychography
%  requires much more signal (imaging dose) than conventional ptychography 
% But the final tomography has still the same number of voxels -> if we search 
% for the voxel values in the tomogram instead of the projections, there
% should be no loss in the signal to noise ratio -> resolution identical to
% the single layer ptychography can be reached. 
% Additionally, if the illumination angle (NA) or the collected NA is larger than
% the tomographic angular sampling angle according to the Crowther criterion 
% It is possible to reconstruct angularly undersampled tomogram and fill
% the gaps by information provided by multilayer ptychography. 
% Of course, the total imaging dose has to be preserved, so the
% ptychography  scanning step or exposure would have to be increased
% 
% Conclusion: eDOF ptychography can benefit from the iterative
% ptycho-tomography 
%
% ================
% Implementation 
% ================
% 1) diffraction from a 3D phantom are created using ptychography code,
% using the simulated sample thickness and other parameters 
% 2) Provide initial 2D reconstructions -> use the fact that computing
% initial guess is cheap -> use conventional ptychtomography as initial
% guess before any iterative ptycho-tomo is started 
% Also and maybe even more importantly, the real measured projections have
% to be aligned first 
% 3) Start to refine the initial guess using the ptycho projections 
%
% P0_i = model complex-valued projection at angle 'i'
% P0c_i = model log complex projection, ie P0 = exp(POc)
% P_i = projection after ptychography update 
% Pc_i = log complex projection after ptychography update 
% Vol = 3D tomographic volume 
% step = step size in the gradient descent method 
% 
% Pseudo code for the iterative ptycho tomo 
% ==========================================
% for i = 1:Nangles
%   Vol_rot = rotate(Vol, theta(i))
%   P0c_i = sum(Vol_rot)       % only in case of single layer, see fwd_proj.m
%   P0 = exp(P0c_i)            % see prepare_distributed_data.m
%   P_i = ptychography(P0_i)   % see ptycho_solver_distributed.m
%   Pc_i = log(Pc_i)           % see prepare_projections.m , note that log is calculated from complex number -> some unwrapping is needed 
%   Vol_upd =  repmat(P_i - P_0, N)  % only in case of single layer ptycho, see back_proj.m 
%   Vol = Vol  + step * rotate(Vol_upd, -theta(i))  % see update_volume.m 
% end
% 
% In case of multilayer, the sum and repmat are replace the centered
% Fourier interpolation. This allows to optimally up/down sample the
% information from different sample layers into a few ptychography layers 
% Therefore, there is no need of propagation between all layers  (ie. size(volData,1))
% of the sample.
%
%
% Important notes: 
%     a) ptychography should run only single iteration and return results back
%     to tomography. But the overhead would be huge -> running a few iterations
%     is better deal. But too many iteration will lead to deterioration of the
%     achievable quality 
%     
%     b) the volume is kept complex values but logarithmized, ie. projection
%     through the volume is exp(sum(volume)) or prod(exp(volume))
%     This is critical for rotation of the volume where some interpolation needs
%     to be used, for example bilinear interpolation. Without the logarithmized
%     values, only nearest neighbours interpolation would be valid 
%     
%     c) Convergence of the gradient descent scheme is slow. Best results are
%     achieved when just a small groups of close to orthogonal angles are applied 
%     in parallel like the SART tomo method. Applying all in parallel like SIRT 
%     would be too slow. 
%     In order to further accelerate the convergence, Nesterov acceleration method 
%     is applied in each iteration to update the tomography volume 
%     
%     d) major issues are GPU memory limitations. Ptychotomography needs to keep 
%     entire tomographic volume in memory, otherwise the data transfer overhead 
%     would significanlty reduce the performance. This limits the maximal size of
%     the volume and some compression schemes, ie store the volume in uint8 or uint16 
%     will be required for datasets large than ~10GB
%     
%     e) ptychography solvers needs to run asynchronously and in parallel with the 
%     tomography solver for maximal performance. Also they have to run on different 
%     GPUs to avoid memory collisions. This asynchronity introduces some noise into the 
%     reconstruction convergence, so the size of the blocks should be kept small. 
%     On the other hand, larger groups uses the resources more efficient. 
% 
%     



%% set shared parameters for all test scripts 
cd(fullfile(fileparts(mfilename('fullpath')), '../'))


addpath('./')
addpath(core.find_base_package)
addpath('./utils')
addpath('../cSAXS_matlab_tomo/')
addpath('../cSAXS_matlab_tomo/utils/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% user settings %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
undersampling = 4; 
Nangles = ceil(300*pi/2 / max(1,undersampling));  % number tomo angles
asize = [192 192]; 
base_path = '../'; 
GPU_id = 1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if undersampling > 1
    warning('Tomogram may be angularly undersampled')
end

%% run common template for initalization
ptycho_path = fileparts(mfilename('fullpath')); 
ptycho_path = replace(ptycho_path, 'tests', ''); 
run(fullfile( ptycho_path, 'template_ptycho.m'))
if ~exist('temporal_data_path', 'var'); temporal_data_path = './temp/'; end
gpuDevice(GPU_id); 
utils.verbose(0)
utils.verbose(struct('prefix', 'ptychotomo'))
utils.report_GPU_usage()

%% general settings
p = rmfield(p, 'engines'); 
p.   scan_number = [];
p.   asize = asize; 
p.   verbose_level = 2;
p.   src_metadata = 'artificial';                                 % load meta data from file; currently only 'spec' is supported;
p.   artificial_data_file = 'tests/test_ML3D_data.m';     % artificial data parameters 
p.   data_prep = 'matlab';                                  % data preparator; 'python' or 'matlab' 
p.   base_path = base_path;
p.   gpu_id = GPU_id; 
%% plotting
p.   save.external = false;                           % Use a new Matlab session to run save final figures (saves ~6s per reconstruction). Please be aware that this might lead to an accumulation of Matlab sessions if your single reconstruction is very fast.
p.   plot.prepared_data = false;                         % plot prepared data
p.   save.store_images = 1;                                       % Write nice jpegs in [p.base_path,'analysis/online/ptycho/'] if p.use_display = 0 then the figures are opened invisible in order to create the nice layout. It writes images in analysis/online/ptycho
p.   plot.fov_box = 1;                                      % Plot the scanning FOV box on the object (both phase and amplitude)
p.   plot.log_scal = [1 1];                                      % Plot on log scale for x and y
p.   plot.positions = 1;                                    % Plot the scanning positions
p.   plot.interval = 50; 

p.   io.data_compression = 0;                                   % file compression for HDF5 files; 0 for no compression
p.   io.file_compression = 0;                                   % file compression for HDF5 files; 0 for no compression


%% SHARE MULTIPLE OBJECTS 
p. share_object = 0;          % Share object between scans. Can be either a number/boolean or a list of numbers, specifying the object index; e.g. [1 2 2] to share the objects between the second and third scan. 


run(fullfile( ptycho_path,p.artificial_data_file))


%% io
p.   ptycho_matlab_path = '';                               % cSAXS ptycho package path
p.   cSAXS_matlab_path = fullfile(ptycho_path, '../cSAXS_matlab_base');                               % cSAXS package path
p.   prepare_data_path = temporal_data_path; 
p.   use_display  = false; 


%% ENGINES
% --------- GPU engines  ------------- 

eng = struct();
eng. name = 'GPU';    
eng. use_gpu = true;            % if false, run CPU code, but it will get very slow 
eng. keep_on_gpu = true;        % keep data + projections on GPU, false is useful for large data if DM is used
eng. compress_data = true;      % use automatic online memory compression to limit meed of GPU memory
eng. gpu_id = [];               % default GPU id, [] means choosen by matlab
eng. check_gpu_load = true;     % check available GPU memory before starting GPU engines 

%% general 
eng. number_iterations = 100;   % number of iterations for selected method 
eng. downscale = 1;             % Ntimes downsize data to make low res. fast guess, similar to presolver engine 
%eng. share_probe = 1;           % Share probe between scans. Can be either a number/boolean or a list of numbers, specifying the probe index; e.g. [1 2 2] to share the probes between the second and third scan. 
eng. method = 'MLs';            % choose GPU solver: DM, RAAR, ePIE, pPIE, hPIE, MLc, Mls, -- recommended are MLc and MLs
eng. opt_errmetric = 'L1' ;     % optimization likelihood - poisson, L1
eng. grouping = 200;            % size of processed blocks, larger blocks need more memory but they use GPU more effeciently
                                % for hPIE, ePIE, MLs methods smaller blocks lead to faster convergence, 
                                % for pPIE, MLc the convergence is similar 
                                % for DM, RAAR is has no effect on convergence
%eng. probe_modes  = 1;          % Number of coherent modes for probe
eng. object_change_start = 1;    % Start updating object at this iteration number
eng. probe_change_start = inf;     % Start updating probe at this iteration number


% regularizations
eng. reg_mu = 0;                % Regularization constant ( = 0 for no regularization)
eng. background = 1e-3; 
eng. delta = 0;                 % press values to zero out of the illumination area, usually 1e-2 is enough 
eng. positivity_constraint_object = 0; % enforce weak positivity in object, usually 1e-2 is already enough 
eng. regularize_layers = 0.01;         % 0<R<<1 -> apply regularization on the reconstructed object layers, 0 == no regularization 

eng. apply_multimodal_update = false; % apply all incoherent modes to object, it can cause isses if the modes collect some crap 
eng. probe_backpropagate = 0;         % backpropagate the probe mask, inf == farfield 

% basic recontruction parameters 
% PIE / ML methods
eng. beta_object = 1;           % object step size, larger == faster convergence, smaller == more robust, should not exceed 1
eng. beta_probe = 1;            % probe step size, larger == faster convergence, smaller == more robust, should not exceed 1
eng. delta_p = 0.1;             % LSQ dumping constant, 0 == no preconditioner, 0.1 is usually safe, 
% DM
eng. pfft_relaxation = 0.1;     % Relaxation in the Fourier domain projection, = 0  for full projection 
eng. probe_regularization = 0.1;% Weigth factor for the probe update (inertia)


% ADVANCED OPTIONS   
% position refinement 
eng. apply_subpix_shift = false;       % apply FFT-based subpixel shift, important for good position refinement but it is slow
eng. probe_pos_search = inf;           % reconstruct probe positions, from iteration == probe_pos_search, assume they are independed
eng. probe_geometry_search = inf;      % reconstruct probe positions, from iteration == probe_geometry_search, assume they have to match geometry model with error less than probe_position_error_max
eng. probe_position_error_max = 10e-9; % max expected random position error of the stages 

% wavefront refinement                  
eng. probe_fourier_shift_search = inf; % refine farfield position of the beam (ie angle) from iteration == probe_fourier_shift_search
eng. estimate_NF_distance = inf;       % try to estimate the nearfield propagation distance  
eng. variable_probe = false;           % Use SVD to account for variable illumination during a single (coupled) scan
eng. variable_SVD_modes = 3;           % OPRP settings , number of SVD modes, apply only for PIE methods 
eng. variable_probe_smooth = 0;        % OPRP settings , apply assumption of smooth evolution of the OPRP modes (ie slow drifts)
eng. variable_intensity = false;       % account to changes in probe intensity


eng_0 = eng; 

addpath(p.   cSAXS_matlab_path)



if 1
    eng = eng_0; 
    eng. method = 'DM';            % choose GPU solver: DM, RAAR, ePIE, pPIE, hPIE, MLc, Mls, -- recommended are MLc and MLs
    eng. grouping = inf; 
    eng. probe_support_radius  = inf; 
	eng. number_iterations = 50;
    eng. probe_change_start = inf; 

    eng. use_display  = true;
    eng. probe_change_start = 10; 
    eng. plot_results_every = inf; 
    [p, ~] = core.append_engine(p, eng);    % Adds this engine to the reconstruction process
    
    eng. method = 'MLc';            % choose GPU solver: DM, RAAR, ePIE, pPIE, hPIE, MLc, Mls, -- recommended are MLc and MLs
    eng. grouping = inf; 
    eng. probe_support_radius  = inf; 
	eng. number_iterations = 50;
    eng. accelerated_gradients_start = 2; 
    eng. use_display  = true;
    eng. probe_change_start = inf; 
    eng. plot_results_every = inf; 
    [p, ~] = core.append_engine(p, eng);    % Adds this engine to the reconstruction process
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE DATA AND INITIAL GUESS %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.prepare.force_preparation_data = true; 
p.use_display = 0; 
angles = 1+ linspace(0,360-1/Nangles, Nangles); 
utils.verbose(-1)
disp('GENERATE DATA')
objects_prepared = {};
for ii = 1:Nangles
   % reset(gpuDevice)
    utils.progressbar(ii, Nangles, 30)
    rng(ii)
    p_tmp = p; 
    p_tmp.   use_gpu = true;
    p_tmp.   verbose_level = -1;
    p_tmp.   positions_pad = round(rand(1,2)*10); 
    p_tmp.   scan_number = ii;
    p_tmp.   rotation_angle = angles(ii); 
    p_tmp.   prepare_data_path = fullfile(ptycho_path,'temp'); 
    
    [pout] = core.initialize_ptycho(p_tmp);
    
    
    core.prep_h5data(pout);
    position_offset(ii,:) =  pout.positions_pad; 
    assert(~core.check_prepared_data( pout ), 'Generated data checked')


    
    objects_ideal{ii} = prod(pout.simulation.obj{1},4);

    p_tmp = pout; 
    % get all engines

    utils.verbose(struct('prefix', {'ptycho'}))
    for ieng=1:length(p_tmp.engines)
        % engine call
        utils.verbose(1, 'Calling engine %s', p_tmp.engines{ieng}.name)
        utils.verbose(struct('prefix', {p_tmp.engines{ieng}.name}))
        [p_tmp, fdb] = core.run_engine(p_tmp,ieng);

        utils.verbose(struct('prefix', {'ptycho'}))
    end
    
    



   objects_prepared{ii} = struct('object', p_tmp.object{1},...
                                 'probe', p_tmp.probes, ...
                                 'positions', p_tmp.positions, ...
                                 'illum_sum', p_tmp.illum_sum{1},  ...
                                 'angle', angles(ii), ...
                                 'position_offset', position_offset(ii,:)); 

    %% plot current projection 
    objects{ii} = p_tmp.object{1};
    o = utils.crop_pad(utils.imshift_fft(p_tmp.object{1}, p_tmp.positions_pad) ,[300,500]); 
    subplot(1,2,1)
    plotting.imagesc3D( abs( o ))
    title('Amplitude')
    colorbar 
    caxis([0.4,1.2])
    colormap bone 
    axis image 
    subplot(1,2,2)
    plotting.imagesc3D( angle( o ))
    colorbar 
    caxis([-1,2])
    title('Phase')
    colormap bone 
    axis image 
    plotting.suptitle(sprintf('Showing initial guess for projection %i', ii))
    drawnow 


end
Npix_simulated_obj = size(pout.simulation.obj{1}); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  TEST FILTERED BACKPROPAGATION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 addpath('../cSAXS_matlab_tomo/')

 obj_size_max = [0, 0]; 
for ii = 1:Nangles 
    obj_size_max = max(obj_size_max, size(objects_prepared{ii}.object)); 
end
     
sino = zeros([obj_size_max,Nangles], 'single'); 
theta = zeros(Nangles,1); 
for ii = 1:Nangles
    utils.progressbar(ii, Nangles);
    sino(:,:,ii) = utils.crop_pad(utils.imshift_fft(objects_ideal{ii},objects_prepared{ii}.position_offset), obj_size_max); 
    theta(ii) = objects_prepared{ii}.angle;
end
% 
sino = utils.crop_pad(sino, floor((obj_size_max - p.asize)/2)*2); 

sino = sino(end/2 + [-20:20],:,:);  % take only a few layers 

par_tomo.pixel_size = pout.dx_spec(1); 
par_tomo.lambda = pout.lambda; 
par_tomo.GPU_list = 1; 
par_tomo.thickness =  p_tmp.simulation.thickness  ; 

sino = gpuArray(sino);
[rec] = tomo.FBP_propagation(sino, theta, 'phase', par_tomo, 0, par_tomo.thickness);

plotting.smart_figure(2322)
% plotting.imagesc3D(rec, 'init_frame', ceil(size(rec,3)/2));
plotting.imagesc3D(max(0,rec), 'init_frame', ceil(size(rec,3)/2));
colormap bone
caxis([0,0.1])
title('Filtered backpropagation reconstruction')
axis image 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test iterative ptychotomo %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


p0 = p;
p0.object_size = pout.object_size; 
p0.asize = p.asize;
p0.simulation.thickness = pout.simulation.thickness;
p0.probes = pout.probes; 




% Make initial object guess 
disp('Store ideal reconstructions')
volData_c = single(squeeze(pout.simulation.obj{1})); 
volData_c = utils.imshift_fft(volData_c,pout.positions_pad);
volData_c = permute(volData_c,[3,2,1]);
volData_c = rot90(volData_c,-1);  % get the same orientation as seen by this code ; 
volData_c = complex(log(abs(volData_c)) , angle(volData_c)); % store in log space instead of complex space 
%%%%%%%%%
Npx_vol = size(volData_c)-100;  % keep only the relevant volume 
volData_c = utils.crop_pad(volData_c, max(Npx_vol(1:2))*ones(1,2));
%%%%%%%%%

% rotate back to the initial angle 
volData_c = gather(utils.imrotate_ax(gpuArray(volData_c), -p_tmp.rotation_angle, 3)); 


plotting.smart_figure(232)
plotting.imagesc3D(max(0,-imag(volData_c)), 'init_frame', ceil(size(volData_c,3)/2));
colormap bone
colorbar
axis image
title('Ideal reconstruction')


%% unwrap data 

addpath('../cSAXS_matlab_tomo/')
objects = objects(1:Nangles);
position_offset = position_offset(1:Nangles, :);
obj_size_max = [0, 0]; 
for ii = 1:Nangles 
    obj_size_max = max(obj_size_max, size(objects{ii})); 
end
for ii = 1:Nangles 
    objects{ii} = utils.crop_pad(objects{ii}, obj_size_max); 
end


disp('Unwrapping ... ')
Npx_proj = [obj_size_max,1];
for ii = 1:Nangles
    utils.progressbar(ii, Nangles,30);
    objects_prepared{ii}.object = utils.stabilize_phase( objects{ii}, 'weights', abs(objects{ii}).^2 ,'fourier_guess', false);
    objects_prepared{ii}.illum_sum = utils.crop_pad(objects_prepared{ii}.illum_sum, Npx_proj);
    delta = 0.01*max(objects_prepared{ii}.illum_sum(:));
    W = sqrt(objects_prepared{ii}.illum_sum.^2 ./ (objects_prepared{ii}.illum_sum.^2 + delta^2)); 
    [object_c] = ptychotomo.prepare_projections(objects_prepared{ii}.object, Npx_proj, p0.asize,false,[], W);
    
    objects_prepared{ii}.object_c = gather(object_c);
    objects_prepared{ii}.weight = W;
    
%     plotting.smart_figure(2454)
%     subplot(1,2,1)
%     plotting.imagesc3D(exp(object_c)); axis off image 
%     title('exp(object_c) - linearized reconstruction', 'interpreter', 'none')
%     subplot(1,2,2)
%     plotting.imagesc3D(objects_prepared{ii}.object); axis off image 
%     title('Original reconstruction')
%     plotting.suptitle('Compare complex projections after unwrapping')
%     drawnow 
    

end
clear object_c


%%

disp('Removing phase offset from initial projection guess')
volData_c = gpuArray(volData_c);
probe_size = p.asize; 
win = tukeywin(Npx_proj(2), probe_size(2)/Npx_proj(2)/2)'.*tukeywin(Npx_proj(1), probe_size(1)/Npx_proj(1)/2);

for ii = 1:Nangles
    utils.progressbar(ii, Nangles);


    [~,objects_prepared{ii}] = ptychotomo.prepare_distributed_data(pout, volData_c, objects_prepared{ii}, [], struct(),  -1, false);
    objects_prepared{ii}.object_c = gather(objects_prepared{ii}.object_c);
    
%     plotting.smart_figure(2455)
%     subplot(1,2,1)
%     plotting.imagesc3D(exp( objects_prepared{ii}.object_c)); axis off image 
%     subplot(1,2,2)
%     plotting.imagesc3D(objects_prepared{ii}.object); axis off image 
%     drawnow 


end

scan_ids = 1:Nangles;
for ii = 1:Nangles
    objects_prepared{ii}.scan_id = scan_ids(ii); 
    objects_prepared{ii}.proj_id = ii; 
    % shift by 1 pixel is needed get the right center of rotation in the projection 
    objects_prepared{ii}.position_offset = objects_prepared{ii}.position_offset - 1; 
end


volData_c = gather(volData_c); 



pout.fmag = []; 
pout.fmask = []; 
pout.simulation.obj = [];
volData0_c = volData_c; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN ITERATIVE PTYCHO-TOMO RECONSTRUCTION  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


par.queue_path = 'reconstruction'; 

% reconstruct the data 
addpath('../cSAXS_matlab_tomo/')
addpath('../cSAXS_matlab_base/')

% smooth volData_c -> make a poor initial guess 
volData_c = utils.imgaussfilt3_fft(volData0_c, 10);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  BASIC SETTINGS %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
par.debug = false; 
par.downsample_angles = 1; 
par.Niter = 10;
par.Niter_inner = 1; 

%% set starting iteration for different actions 
par.start_3D_reconstruction = 1;  % start reconstructing the 3D volume 
par.ptycho_interval =  1;  % solver ptycho each N-th call 
par.ptycho_reconstruct_start = 2;   % start ptychography reconstructiom 
par.ptycho_ML_reconstruct_start = par.ptycho_reconstruct_start; %  start multilayer solver 
par.ptycho_accel_start =  par.ptycho_ML_reconstruct_start+1 ;  % start momemntum acceleration 
par.plot_every = 10;   % seconds


par.Nlayers_max = 32; % maximal number of reconstructed layers 
par.apply_support = true; % apply support constraint around the reconstructed volume 
par.smooth_reconstruction = false;  % apply smoothing filter on the reconstruction, helps agains artefacts fom the rotation 
par.wait_time_solver = 10; % how long should the code wait for a projection before giving up 
par.max_queue_length = 5;  % how many projection should be kept in the processing queue. More == better parallelism but more unstable solver 

par.lambda = 0.5;  % relative update step length in gradient descent 

par.prepare_data_path = './temp/S%05i/'; 


%% set constraints for the 3D reconstruction, it makes solver more stable 
par.max_value = 0 + 0i; 
par.min_value = -4.0e-3 -10e-2i;  
 

%%%%%%%%%%%%%%%%%%%%%%%%%
%  set ptychography     %
%%%%%%%%%%%%%%%%%%%%%%%%%

p0.verbose_level = 0;

% set data preparation 
p0.queue.name = 'filelist'; 
p0.queue.file_queue_timeout = 0.1; 
p0.queue.recon_latest_first = false; 

p0.prepare.auto_prepare_data = true; 
p0.prepare.force_preparation_data = false;

% set reconstruction engines 
p0.number_iterations = 5;
p0.plot.interval = inf; 
p0.probe_change_start = 1;
p0.object_change_start = 2;


% engine settings 
p0.engines = {struct()};
p0.engines{1}. name = 'GPU';   
p0.engines{1}. method = 'MLc';   

p0.engines{1}.verbose_level = 0;
p0.engines{1}.mirror_objects = false; 
p0.engines{1}.share_object = false;
p0.engines{1}.probe_support_fft = true;            % Apply probe support in Fourier space, ! uses model zoneplate settings to estimate support size 
p0.engines{1}.probe_support_radius  = []; 
p0.engines{1}.grouping = inf;   % group size 
p0.engines{1}.regularize_layers = 0;
p0.engines{1}.use_display  = false;
p0.engines{1}.momentum = 0;
p0.engines{1}.probe_support_fft = false;
p0.engines{1}.delta = 0;
p0.engines{1}.beta_LSQ = 0.9;
p0.engines{1}.delta_p = 0.1;
p0.engines{1}.remove_object_ambiguity = false;
p0.engines{1}.accelerated_gradients_start = 2;

p0.prepare_data_filename = '';
p0.scan = '';
p0.fmag = []; 
p0.fmask = []; 
p0.simulation.obj = [];
p0.object_size = size(objects_prepared{1}.object); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p0.engines{1}.probe_pos_search = inf;





Npx_vol = size(volData_c); 

volData_c = utils.crop_pad(volData_c, ceil(Npx_vol / 64) * 64); 

% correct for updated volume size -> keep thicnkness per pixel constant 
p0.thickness = pout.simulation.thickness * size(volData_c,1) /  Npix_simulated_obj(4); 

utils.verbose(0)
utils.verbose(0, 'p0.simulation.thickness %g microns', p0.simulation.thickness*1e6)




par.force_initialization = true; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% START PTYCHO SOLVERS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run asynchronous solver in background , there can be many of them
% started, they will automatically use different GPUs


Nsolver = 1;   % how many solvers should be started 
for ii = 1:Nsolver
    system('matlab -nodisplay -r ptychotomo.call_tomo_solver_asynchronous & '); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% START TOMO SOLVER  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% run tomo solver in foreground 
volData_rec = ptychotomo.tomo_solver_distributed(par, volData_c, objects_prepared, p0, angles,scan_ids); 


% clear temporal data 
rmdir('./temp/', 's')
rmdir('./reconstruction/', 's')

