function [p_out, fdb] = GPU_MS(p)

    import math.*
    import utils.*
    import plotting.*
    import engines.GPU_MS.*
        
    verbose(struct('prefix','GPU/CPU_MS-engine'))

    %%%%%%%%%%%%%%%%%%  WRAPPER FOR FOR GPU CODE %%%%%%%%%%%%%%%%%%%

    if p.verbose_level > 2 && check_option(p, 'meta') && isfield(p.meta{1}, 'spec')
        verbose(1, 'Starting GPU_MS engines on scans:')
        for i = 1:length(p.meta)
            verbose(1,p.meta{i}.spec.S)
        end
    end    
       
    fdb = []; 

    %% load default settings 
    param = initialize.get_defaults();
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% load values from p-struct to self-class and param structure %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [self, param, p] = initialize.load_from_p(param, p);
    
    global use_gpu gpu 
    verbose_level = param.verbose_level;
  
    %% INITIALIZE GPU IF AVAILABLE %%% 
    gpu_id = -1; 
    param= GPU_wrapper.initialize(param);
    use_gpu = param.use_gpu;

    
    if use_gpu 
        gpu_id = gpu.Index;
        verbose(struct('prefix','GPU'))
    else
        verbose(struct('prefix','CPU'))
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%%%%%%%%% GPU SOLVER %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    time_start = tic; 
    
    try
        
        %% initialization of the matlab engine 
        
        % rescale p -> for low resolution initial guess 
        tmp = initialize.rescale_inputs(self, param.Np_p_presolve, true);
        if param.mirror_objects
            tmp = shared.flip_mirror_scan(tmp, true); 
            param.share_probe = false;  % probe will be flipped as well !!
        elseif param.align_shared_objects && size(tmp.object,1) > 1 && param.share_object
            tmp = shared.align_objects(tmp); 
        end
               
        % make final preparations & checks before engine is loaded 
        [tmp,param] = initialize.check_inputs(tmp, param);   
        
        %disp(size(tmp.probe))

        [tmp,cache] = initialize.init_solver(tmp, param);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% WAIT IN QUEUE UNTIL THERE IS ENOUGH FREE GPU MEMORY TO RUN
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Check RAM 
        [mem_avail, mem_total] = utils.check_available_memory(); 
        verbose(1,'%.0f of %.0f GB RAM free', mem_avail/1e3, mem_total/1e3 ); 

        while  param.use_gpu
            if isinf(param.grouping)
                grouping_tmp = 5;  % for inf grouping, use a small initial guess of the group size and refine later for given availible memory 
            else
                grouping_tmp = param.grouping; 
            end
            [required_mem , data_mem, object_mem, fft_mem] =  GPU_wrapper.estimate_req_memory(tmp, param, grouping_tmp); 


            % get an estimate of the free memory 
            if param.check_gpu_load
                Ntest = 5; % [s]
                verbose(1,'Check GPU load')
            else
                Ntest = 1; 
            end
            if utils.verbose >= 0 && gpu.TotalMemory - gpu.AvailableMemory > 1e9
                % if someone uses the GPU, report the user
                %utils.report_GPU_usage %disabled by YJ to avoid potential error
            end

            for ii = 1:Ntest
                if param.check_gpu_load && verbose_level > 1; progressbar(ii, Ntest); end
                available_mem(ii) = gpu.AvailableMemory;
                if isnan(available_mem(ii)); verbose(0, 'GPU reset'); reset(gpu); continue; end
                if available_mem(ii) > 2*required_mem; break; end 
                pause(1)
            end

            if param.check_gpu_load && verbose_level > 1
                progressbar(Ntest, Ntest)
                fprintf('\b'); 
            end
            if ~isinf(param.grouping)
                verbose(1,'\n----MEMORY REPORT GPU ID:%i-----\nRequired total memory %3.2fGB\n - Required data memory %3.2fGB \n - Required object memory %3.2fGB \n - Required FFT memory %3.2fGB \n ============================ \navailable memory %3.2f/%3.2fGB',...
                    gpu_id, required_mem / 1e9, data_mem / 1e9, object_mem / 1e9, fft_mem/1e9, min(available_mem)/ 1e9, gpu.TotalMemory/ 1e9)
            end
            if required_mem < min(available_mem)
                % ready to go .... 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%% call GPU solver  %%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                verbose(1,' === Starting %s solver === ', param.method)

                try
                    % call the solver          
                     [out, fourier_error, fdb.score] = ptycho_solver(tmp, param, cache);
                    
                catch ME
                    if any(strcmpi(ME.identifier,{'id:parallel:gpu:array:OOMForOperation','id:MATLAB:LowGPUMem','MATLAB:LowGPUMem',...
                            'parallel:gpu:array:OOM','parallel:gpu:device:UnknownCUDAError', ...
                            'parallel:gpu:array:OOMForOperation', 'parallel:gpu:device:UnknownCUDAError', ...
                            'parallel:gpu:array:FFTInternalError'}))
                        warning('\n=== Failed due to GPU issue trying again ... === \nFile %s Line: %i \n  id:%s msg:%s \n \nGPU ID:%i\nRequired memory %3.2fGB\navailable memory %3.2f/%3.2fGB\n', ...
                        ME.stack(1).name,  ME.stack(1).line, ME.identifier, ME.message, gpu_id, required_mem / 1e9, gpu.AvailableMemory/ 1e9, gpu.TotalMemory/ 1e9' )
                        % try to return back to waiting queue
                        
                        verbose(-1,'Reset GPU')
                        try
                            gpuDevice(gpu_id)
                        end
                        verbose(-1,'Return back to queue')
                        required_mem = required_mem .* 1.2; % assume that the required memory was too low, -> increase 
%                         keyboard;
                        continue
                    else
                        rethrow(ME)
                    end
                end
                break
            elseif required_mem > gpu.TotalMemory * 0.9
                error('Too large memory requirements for selected GPU, try to reduce grouping')
            end
            % otherwise keep waiting 
            wait_time = 5;
            warning('Low memory on GPU %i, waiting %is ...',gpu_id, wait_time)
            pause(wait_time)
            
        end
        if ~param.use_gpu
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%% call CPU solver (backup) %%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [out, fourier_error, fdb.score] = ptycho_solver(tmp, param, cache);
        end
        
    catch ME
        try
            warning('\nFile %s Line: %i id: %s msg:%s \n \nGPU ID:%i\nRequired memory %3.2fGB\navailable memory %3.2f/%3.2fGB\n', ...
                ME.stack(1).name,  ME.stack(1).line,ME.identifier, ME.message, gpu_id, required_mem / 1e9, gpu.AvailableMemory/ 1e9, gpu.TotalMemory/ 1e9' )
        catch
            verbose(0,'File %s Line: %i \n id:%s msg:%s ', ME.stack(1).name,  ME.stack(1).line, ME.identifier, ME.message)
 
        end
        if verbose < 1
            rethrow(ME)
        else
            fprintf('verbose %i', verbose)
           keyboard 
        end
    end    
    
    if check_option(p, 'clean_residua')
        warning('Cleaning residua, it may introduce errors in lowest spatial frequencies')
        out = engines.GPU.shared.clean_residua(out, cache); 
    end
    if param.mirror_objects
        % flip the object, data and positions back 
        out = shared.flip_mirror_scan(out, false); 
    end
    %% upscale outputs  back if needed and update valued in "self" structure

    out.diffraction = self.diffraction; 
    out.mask = self.mask; 

    out = initialize.rescale_inputs(out, p.asize, false);

    % do not return the modified data, return the original 
    p_out = initialize.save_to_p(out, param, p, fourier_error); 

    verbose(0,' === Finished %s solver === in %4.3gs', param.method, toc(time_start))

    fdb.status.status = false;   
    % reset verbosity back to the p-struct level 
    verbose(p.verbose_level)
    

    
end

