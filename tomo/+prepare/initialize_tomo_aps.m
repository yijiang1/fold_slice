% INITIALIZE_TOMO_APS basic initialization steps of tomography -> check validity of the inputs, 
% load first projection and store its parameters, check angles, create output folders 
% Created by YJ Based on PSI's function

%  [par, angles_check, object] = initialize_tomo_aps(par, scans, use_gpu, object_preprocess_fun)
%
% Inputs: 
%   **par         - basic parameters defined in template 
%   **scans       - list of the scans to be loaded 
%   **use_gpu     - (bool), dont use GPU if use_gpu == 0, (default = true )
%   **object_preprocess_fun - user defined preprocessing function applied on the loaded projections, e.g. in laminography it can be rotation, default = @(x)x
%
% *returns* 
%   ++par           updated basic parameters 
%   ++angles_check	 
%   ++object        example of one loaded projection 


function [par, angles_check, object] = initialize_tomo_aps(par, scans, use_gpu, object_preprocess_fun)

    import ptycho.*
    import io.*
    utils.verbose(struct('prefix', 'initialize'))
    %% initial checks 
    if verLessThan('matlab', '9.3')
        warning on
        warning('Only Matlab versions >= 2018a are tested and supported, \nYour Matlab version is %s', version) 
        pause(5)
    end
    if nargin < 3
        use_gpu = true;
    end
        
    if gpuDeviceCount == 0 && use_gpu 
       warning('Using CUDA enabled GPU is strongly recommended') 
       pause(5)
       use_gpu = false;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%   CHECK GPU AVAILIBILITY %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if use_gpu
        if gpuDeviceCount == 0
           error('Code needs CUDA enabled GPU, suppress by setting input parameter "use_gpu=false" ') 
        end
        if any(par.GPU_list > gpuDeviceCount)
           error('Selected GPU in GPU_list is not available')
        end
        gpu = gpuDevice(par.GPU_list(1)); 
        if ~verLessThan('matlab', '9.4') &&  gpu.ToolkitVersion < 9
           error('Code needs CUDA 9.0 to work with Matlab 2018a and newer')
        elseif verLessThan('matlab', '9.4') &&  gpu.ToolkitVersion < 8
           error('Code needs at least CUDA 8.0 to work with Matlab 2017b') 
        end
        fprintf('=================================================== \n')
        fprintf('=== Available memory for GPU %i : %2.1fGB / %2.1fGB === \n', gpu.Index, gpu.AvailableMemory/1e9, gpu.TotalMemory/1e9)
        fprintf('=================================================== \n')
        
        % check that more than 3GB of GPU mem is free and that 90% of total
        % memory is available -> make sure that this template is the only
        % process using the selected GPU 
        reset(gpu)
        if ~debug() &&  (gpu.AvailableMemory < gpu.TotalMemory * par.check_gpu_percentage || gpu.AvailableMemory < 3e9)
        
            utils.verbose(0,'\n\n=============== GPU report ================')
            !nvidia-smi
            warning on 
            warning off backtrace 
            if gpu.AvailableMemory < gpu.TotalMemory * 0.9
                warning(['Memory in GPU  %i (Nvidia id:%i) is probably used by other user,'...
                'try to manunally choose GPU or kill other processes'], gpu.Index,  gpu.Index-1)
            else
                warning(['Memory in GPU  %i (Nvidia id:%i) is less than recommended 3GB,'...
                'try to manunally choose GPU or kill other processes'], gpu.Index,  gpu.Index-1)
            end
            warning on 
            
            % check who is using the GPU 
            utils.report_GPU_usage(gpu.Index); 
            
            if ~debug() && ~par.online_tomo
                if  ~strcmpi(input('Do you want to continue [y/N]', 's'), 'y')
                   error('Set other GPU to use by par.GPU_list parameter') 
                end
            end
            % this is only recommende value, the code should run even with
            % less, but then it gets less efficient.           
        elseif (gpu.AvailableMemory < gpu.TotalMemory * 0.9 || gpu.AvailableMemory < 3e9)
            %utils.report_GPU_usage
            
        end
        
    end
    if nargin < 4
        object_preprocess_fun = []; % no preprocessing function 
    end

    par.use_GPU = use_gpu;  % store user preferences in using GPU 
    %{
    %% First check if reconstructions exist 
    proj_file_names = {};

    hasRecon = zeros(length(scans),1);
    for ii=1:length(scans)
        progressbar(ii, length(scans))
    	file = find_ML_recon_files_names(par, scans(ii)); %find ML recon outputs
        if ~isempty(file)
            hasRecon(ii) = 1;
        end        
    end
    %}
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% initial values - LOAD ONE FRAME FOR DEFINING PTYCHO SCAN VALUES %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    file = [];
    ii = 1;
    while isempty(file) && ii <= length(scans)
        %file = find_projection_files_names_aps(par, scans(ii));
        file = find_ML_recon_files_names(par, scans(ii)); %find ML recon outputs
        if isempty(file)
            %warning(['Out of luck - Reconstruction not found']);
            disp(['Out of luck - Reconstruction not found']);
            ii = ii+1;
        else
            break
        end
    end
    if isempty(file)
        error('No reconstructions found, check that analysis folder path contains scans %i-%i', min(scans), max(scans))
    end
    disp(file)
    %%%  Read first projection to check size and reconstruction parameters
    display(['Reading file: ' file])
    [object, probe, dx_spec] = load_aps_ML_recons(file);
    if iscell(probe)
        probe = single(probe{1}(:,:,1));   % keep only the first mode
    else
        probe = single(probe(:,:,1));   % keep only the first mode
        
    end
    par.asize = size(probe); % probe size

    par.dims_ob_loaded = [size(object,1), size(object,2)];  % load the sizes directly from the object, note that "object_preprocess_fun" can crop/rotate the image !!
    %{
    if isfield(p, 'scanindexrange')
        p.scanidxs{1} = p.scanindexrange(1):p.scanindexrange(2);
        positions = int32(p.positions(p.scanidxs{1},:));
        indices = int32(1:length(p.scanidxs{1})); 
        % get at least some estimation of the illumination intensity for different regions in the
        % projection 
        par.illum_sum = utils.add_to_3D_projection(abs(probe).^2,zeros(max(p.object_size,[],1),'single'),positions,indices, true); 
    else
        % if nto availible, get et least a crude guess 
        par.illum_sum = ones(par.dims_ob_loaded-par.asize);
    end
    %}
    par.illum_sum = ones(par.dims_ob_loaded-par.asize);
    par.illum_sum = utils.crop_pad(par.illum_sum,par.dims_ob_loaded);
    par.illum_sum = par.illum_sum ./ quantile(par.illum_sum(:), 0.9);  % normalize the values to keep maximum around 1

    % in case of unequal pixel size 
    if dx_spec(1) ~= dx_spec(2)
        % upsample the data in the dimennsion with lower resolution (-> at least relax issues in tomography interpolation)
        pixel_scale = dx_spec ./ min(dx_spec) ; 
        dims_ob_new = round(par.dims_ob_loaded .* pixel_scale); 
        par.illum_sum = max(0,real(utils.interpolateFT(par.illum_sum, dims_ob_new)));
        object = utils.interpolateFT(par.illum_sum, dims_ob_new);
        par.asize = round(par.asize .* pixel_scale);
        probe = utils.interpolateFT(probe, par.asize);
        dx_spec(:) = min(dx_spec);
    end
        
    if ~isempty(object_preprocess_fun)
        % apply custom preprocessing, e.g. rotation and flipping for
        % laminography setup
        object = object_preprocess_fun(object); 
        par.illum_sum = max(0, object_preprocess_fun(par.illum_sum)); 
    end
    
    par.dims_ob = [size(object,1), size(object,2)];  % object size after object_preprocess_fun
    par.probe = probe; 
    if ~isfield(par, 'lambda')
    	par.lambda  = p.lambda;                         % wavelength [m]
    end
    par.pixel_size=dx_spec(1) *  2^par.downsample_projections;       % reconstructed pixel size [m]
    if dx_spec(1)~=dx_spec(2)
        warning('Pixel size not symmetric - This code cannot handle')
    end

    par.factor=par.lambda/(2*pi*par.pixel_size); 
    par.factor_edensity = 1e-30*2*pi/(par.lambda^2*2.81794e-15);
    
    %{
    %%%  Check angles %%%
    if par.checkangles
        [par.scans_check, angles_check] = tomo_angles(projections, subtomograms, ...
            scan_num, subs_to_do); % ignores the repeated 180deg scan.
    else
        angles_check = [];
    end
    %}
    angles_check = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% GENERATE SCAN STRING FOR FILES DESCRIPTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    par.scans_string = {}; 
    if isfield(par, 'output_folder_prefix') && ~isempty(par.output_folder_prefix)
        par.scans_string{end+1} = par.output_folder_prefix; 
    end
    
    if ~isempty(par.tomo_id)
        auxstr = repmat('%i+',1,length(par.tomo_id)); 
        par.scans_string{end+1} = sprintf(['id_',auxstr(1:end-1)], par.tomo_id);
    elseif  par.online_tomo
        par.scans_string{end+1} = sprintf('S%05d',scans(1));
    end
    %{
    % load sample name if provided 
    if ~isfield(p, 'samplename')
        par.samplename = '';
    else
        par.samplename = p.samplename;
    end
    
    if ~isempty(par.samplename)
        par.scans_string{end+1} =  par.samplename; 
    end
    %}
    if ~par.online_tomo
        par.scans_string{end+1}= sprintf('S%05d_to_S%05d',scans(1),scans(end));
    end
    
    par.scans_string = join(par.scans_string, '_');
    par.scans_string = par.scans_string{1}; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Output folder  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    par.output_folder= {fullfile(par.output_path, 'tomo'), par.scans_string, par.filesuffix, par.fileprefix}; 
    if par.online_tomo
        par.output_folder{end+1}= 'online';
    end
    
        
    par.output_folder = join(par.output_folder, '_');
    par.output_folder = par.output_folder{1}; 
    

    %{
    if  ~debug()
        utils.verbose('Output folder: %s', par.output_folder)
        if ~exist(par.output_folder,'dir') 
            mkdir(par.output_folder); 
        end
        [~,attr] = fileattrib(par.output_folder); 
        if ~(attr.UserWrite || attr.GroupWrite)
            error('Output path %s is not writable', par.output_folder)
        end
        % For website
        subdir_online = fullfile(par.base_path,'analysis/online/tomo/');
        if ~exist(subdir_online,'dir')
            mkdir(subdir_online);
        end
        par.online_tomo_path = sprintf('%sonline_tomo_S%05d', subdir_online, min(par.scanstomo));

    end
    %}

end
