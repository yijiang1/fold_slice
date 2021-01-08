%    BLOCK_FUN autosplitting function for fast GPU data processing. 
%    in the simplest case this function is equivalent to 
%           [varargout{:}] = fun(fp16.get(stack_object(ROI{:},:)), varargin{2:end}); 
%    But for large datasets, it will split the stack_object along 3rd axis and try to also
%    split other inputs if possible. Then it will submit each of the
%    blocks on 1 or more GPUs and process. Processed data are returned back to RAM. 
%    it is useful even to process large blocks of data in RAM, to prevent too large RAM allocation. 
%    NOTE: Avoiding GPU processing ('use_gpu', false) can be useful for simple operations when the memory transfer to GPU is much slower than task execution !!
%    NOTE: Multi GPU processing addes extra overhead for moving data to shared arrays and back, useful only for very expensive operations 
%    NOTE: If input is fp16 precision, it is automatically converted to singles for processing and result are converted back to fp16 when returned (unless ('use_fp16', false) is used)
%    NOTE: Automatic splitting assumes memory required for two FFT operation and real to complex conversion with single precision numbers, ie MEMORY = 4*6*8*numel(stack_object)
%    
%    varargout = block_fun(fun, stack_object, varargin)
%
%     Inputs:
%       **fun               - function to be executed 
%       **stack_object      - 3D array that will be split along 3rd dimension and processed on GPU (as default)
%       **varargin          - other arguments, if their size along 1st dimension is the same as 3rd axis of stack_object, they will be split 
%       **parameter structure - if the last argument is structure, 
%               it will not be given to the "fun" but it will be used as a parameter input for block_fun. 
%               SEE DESCRIPTION IN CODE FOR ALL POSSIBLE INPUTS 
%     *returns*
%       ++varargout           - 3D array of size stack_object or 1D/2D with 1st dimension assumed to be equal to 3rd axis of stack_object
%
%     Examples: 
%       % simple 1:1 array processing 
%       X = randn(Nx, Ny, Nz)
%       X_fft = block_fun(@fft2,X);   % equivalent to fft2(x) for small datasets 
%       % simple 1:1 array processing with manually defined blocks size , can be useful if autosplitting does not provide good results (ie too large / too small blocks)
%       X_fft = block_fun(@fft2,X, struct('Nblocks', 10));   % equivalent to fft2(x) for small datasets 
%       % information reduction to (Nz,1,1), 
%       [x,y] = block_fun(@utils.center, X); % equivalent to utils.center(x) for small datasets 
%       % information reduction to (Nx, Nx,1)
%       X_accum = block_fun(@(x)(abs(x).^2), X, struct('reduce_fun', 'plus'));   % equivalent to sum(abs(x).^2,3) for small datasets 
%       % information reduction to (Nx, Nx,1) using only CPU 
%       X_accum = block_fun(@(x)(abs(x).^2), X, struct('use_gpu', false, 'reduce_fun', 'plus'));   % equivalent to sum(abs(x).^2,3) for small datasets but avoids using GPU 
%       % information reduction to (Nx, Nx,1) and return results as single even if inputs is fp16
%       X_accum = block_fun(@(x)(abs(x).^2), X, struct('use_fp16', false, 'reduce_fun', 'plus'));   % equivalent to sum(abs(x).^2,3) for small datasets but avoids using GPU 
%       % use multiple inputs and perform inplace image shifting in order to prevent large memory allocation 
%       shift_x = randn(Nz,1); shift_y = randn(Nz,1); 
%       X_accum = block_fun(@utils.imshift_fft, X, shift_x, shift_y, struct('inplace', true));   % equivalent to imshift_fft(X, shift_x, shift_y)
%       % process only a small subregion of input array, avoid slow matlab memory copy / allocation 
%       ROI = {1:10, 4:50}; 
%       X_phase_cropped = block_fun( @angle, X, struct('ROI', ROI));   % angle(X(ROI{:}))
%       % make the processing quiet, ie no progressbars are shown 
%       X_phase = block_fun( @angle, X, struct('verbose_level', 0));   % angle(X)




% clean all memory blocks:
% ! ipcs -m | cut -d' ' -f2 | grep '^[0-9]' | while read x; do ipcrm -m $x; done


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



function varargout = block_fun(fun, stack_object, varargin)

    global pprev
    assert(~isempty(stack_object), 'Processed array is empty')

    Nlayers = size(stack_object,3); 
    varargout = cell(nargout,1); 
    varargin = [{stack_object}, varargin]; 
    % check if there is given parameter structure
    if isstruct(varargin{end}) 
        param = varargin{end}; 
        varargin = varargin(1:end-1);   % !!! remove the parameter structure from the vargin list
    else
        param = struct();
    end
    ninputs = length(varargin); 
    assert(~isempty(stack_object), 'Input array is empty')
    
    %% ============= LOAD PARAMETERS FROM PARAM STRUCTURE OR USE PREDEFINED DEFAULTS ================
    if ~isfield(param, 'verbose_level')
        param.verbose_level = 1; 
    end
    if ~isfield(param, 'use_GPU')
        param.use_GPU = true;   % some operations are just faster without moving to GPU 
    end  
    if ~isfield(param, 'move_to_GPU')
        param.move_to_GPU = true; % move blocks directly yo GPU or let the called function to decide 
    end 
    if ~isfield(param, 'GPU_list') || isempty(param.GPU_list)
        if param.use_GPU && gpuDeviceCount > 0
            gpu = gpuDevice; 
            param.GPU_list = gpu.Index; 
        else
            param.GPU_list = -1;   % list of the used GPUs
        end
    end 
    if ~isfield(param, 'inplace')
        param.inplace = false;  % run the function inplace , !! dont stop the operation in middle of process !! 
                                % assume that output is the same as first
                                % input array and it will write results
                                % into this memory without reallocating
    end 

    if ~isfield(param, 'ROI') || isempty(param.ROI)   % process only limited ROI of all inputs of size stack_object 
        param.ROI = {':', ':'}; 
    end  
    Nelements = 1;
    for i = 1:2
        if strcmpi(param.ROI{i},':') 
            param.ROI{i} = 1:size(stack_object,i); 
        end
        Nelements = Nelements * length(param.ROI{i});
    end
    assert(Nelements>0, 'Selected ROI is empty')
    
    if ~isfield(param, 'Nblocks')  % number of blocks to split the input volume, [] = auto 
        param.Nblocks = []; 
    end
    if ~isfield(param, 'full_block_size')  % largest volume size in the calculation, default = size(stack_object)
        param.full_block_size = size(stack_object); 
    end
    if ~isfield(param, 'use_fp16')
        param.use_fp16 = isa(stack_object, 'uint16'); 
    end
    % function applied on the computed data to reduce along 3rd axis 
    if ~isfield(param, 'reduce_fun')
        param.reduce_fun = []; 
    else
        assert(ismember(func2str(param.reduce_fun), {'min','max','plus'}))
    end
    if ~isfield(param, 'use_shared_memory')
        param.use_shared_memory = false;  % use shared memory for data exchange, it is used only if N_GPU > 1 or for debugging 
    end
    if param.inplace
        varargout{1} = stack_object; 
    end

    if param.inplace
        warning on
        warning off backtrace
        warning('Inplace data processing, DO NOT INTERRUPT')
        warning on backtrace
    end
    
    if (gpuDeviceCount == 0 || isa(stack_object, 'gpuArray') || ismatrix(stack_object)) && ...
            (isempty(param.Nblocks) || param.Nblocks == 1)
        %% ================ EXECUTE FUN =============================
        [varargout{:}] = fun(stack_object(param.ROI{:},:), varargin{2:end}); 
    else
        %% GPU or CPU splitting is required
        N_GPU = max(1,length(param.GPU_list)); 
        if param.use_GPU
            gpu = gpuDevice; 
            if ~ismember(gpu.Index,param.GPU_list)
                gpu = gpuDevice(param.GPU_list(1)); 
            end
        end
        if isempty(param.Nblocks)
            %% autoestimate block size 
            if param.use_GPU
                %% empirical condition assuming FFT involved, may be too pesimistic 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Nblocks =  ceil( (4*6*8* prod(param.full_block_size)) /  gpu.AvailableMemory) ; 
                Nblocks = max(Nblocks,  prod(param.full_block_size)/ double(intmax('int32'))); 
                Nblocks = max(Nblocks, N_GPU); 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                %% if CPU is used, split it to <20GB  blocks 
                max_block_size = min(utils.check_available_memory*1e6, 20e9); %% work with 10GB blocks 
                Nblocks  = ceil( (6*8* prod(param.full_block_size)) / max_block_size) ; 
            end
        else
            Nblocks = param.Nblocks; 
        end
        
        % optimize size of the block to equally use all provided GPUs 
        N = ceil(Nlayers / Nblocks); 
        Nblocks = ceil(Nlayers / N / N_GPU) * N_GPU;   % make it splitable to N_GPU
        N = ceil(Nlayers / Nblocks);         
               
        % parse inputs and store infomations needed to split the inputs 
        splitable_along_3rd_axis = false(ninputs,1);
        apply_ROI = false(ninputs,1);
        
        for jj = 1:ninputs
            splitable_along_3rd_axis(jj) =  ndims(varargin{jj}) == 3 && (size(varargin{jj},3) == size(stack_object,3)); 
            apply_ROI(jj) = all([size(varargin{jj},1),size(varargin{jj},2)] == [size(stack_object,1), size(stack_object,2)]) && ...
                        ((size(varargin{jj},1) ~= length(param.ROI{1})) || (size(varargin{jj},2) ~= length(param.ROI{2})));
        end
                
        % expected inputs are singles or uint16, if they are doubles, convert to single so save memory and speed up execution on GPU     
        for jj = find(splitable_along_3rd_axis)'
            if ~ismember(class(varargin{jj}), {'single', 'uint16', 'uint8'})
                varargin{jj} = single(varargin{jj}); 
            end
        end
        
        if Nblocks == 1
            %% in case no splitting is needed
            varargout = process_in_single_block(fun, varargin, apply_ROI,splitable_along_3rd_axis, nargout, param);
            return
            
        else
            
            % when the function is finished, make sure to execute following code 
            global status 
            status = true; 

            % avoid using shared memory if not needed, it is a bit slower
            use_shared_memory = length(param.GPU_list) > 1 || param.use_shared_memory;
            if use_shared_memory
                out = onCleanup(@myCleanupFun);
            end
            for ii = 1:length(varargin)
                if splitable_along_3rd_axis(ii) || N_GPU > 1
                   varargin{ii}  = gather(varargin{ii}); 
                end
            end
            
            ind = {};

            for ii = 1:Nblocks
                ind{end+1} = 1+N*(ii-1) : min(Nlayers, N*ii); 
                % avoid single layer arrays 
                if length(ind{end}) < 2
                    ind{end-1} = [ind{end-1}, ind{end}]; 
                    ind(end) = []; 
                end
            end
            Nblocks = length(ind); 

            %% open parpool to allow multi GPU processing 
            if  N_GPU > 1
                poolobj = gcp('nocreate');
                if isempty(poolobj) || poolobj.NumWorkers < N_GPU
                    delete(poolobj);
                    poolobj = parpool(N_GPU);
                end
                poolobj.IdleTimeout = 600; % set idle timeout to 10 hours
            end

            if param.verbose_level>0
                pprev = -1; 
            end
            %% START OF OUTER GPU LOOP 
            outputs_blocks = [];
            %% unitialize one solver per GPU 
            for block_id = 1:N_GPU
                % parse inputs and try to split them if possible 
                [outputs_blocks, inputs_block{block_id}] = submit_block(block_id,block_id, ind, varargin,outputs_blocks ,param, fun, apply_ROI, use_shared_memory, Nlayers,nargout, splitable_along_3rd_axis ); 
            end
            
            unprocessed_blocks = N_GPU+1:Nblocks;
            
            %% merge blocks back from GPUs
            for ii = 1:Nblocks
                if param.verbose_level>0; utils.progressbar(ii, Nblocks); end
                % set values from the small blocks to the final output arrays
                [thread_id, varargout] = gather_block(varargout,param,Nlayers, nargout,ind, N, outputs_blocks); 
                if ~isempty(unprocessed_blocks)
                    block_id = unprocessed_blocks(1);
                    unprocessed_blocks(1) = []; 
                    % submit a new job once the previous is  finished 
                    [outputs_blocks, inputs_block{thread_id}] = submit_block(block_id,thread_id,ind, varargin,outputs_blocks ,param, fun, apply_ROI, use_shared_memory, Nlayers,nargout, splitable_along_3rd_axis ); 
                    if  isa(outputs_blocks, 'parallel.FevalFuture') && any(cat(1,[outputs_blocks.Read]))
                        outputs_blocks
                        warning on 
                        warning('Unknown error in parallel processing')
                        keyboard
                    end
                end

            end

            % END OF OUTER GPU LOOP 
        end
    end

    % everything was fine -> no cleaning needed
    status = false; 

    if param.inplace
        warning on
        warning off backtrace
        warning('Inplace data processing finished')
        warning on backtrace
    end
end

%% auxiliary function 

function [outputs_blocks,inputs_block] = submit_block(block_id, thread_id, ind, inputs, outputs_blocks, param, fun, apply_ROI, use_shared_memory, Nlayers , Noutputs,  splitable_along_3rd_axis )
    inputs_block = inputs;
    ind = ind{block_id}; 
    % get blocks of data 
    for jj = 1:length(inputs)
        if splitable_along_3rd_axis(jj)
            if apply_ROI(jj)
                ROI = param.ROI; 
            else
                ROI = {}; 
            end
            % use fast splitting based on MEX files for large splitable blocks 
            inputs_block{jj} = get_block(inputs_block{jj}, ind, ROI,use_shared_memory); 
        elseif isnumeric(inputs_block{jj}) && size(inputs_block{jj},1) == Nlayers 
            % !! assume split along first dimension of the object !!!
            inputs_block{jj} = inputs_block{jj}(ind,:,:,:);   % split it and provide subfunction only the valid chunk
        end
    end
    if isempty(outputs_blocks); clear outputs_blocks; end
    %% ----------- call the function "fun" , execute on GPUs -----------------
    if length(param.GPU_list) <= 1 || ~param.use_GPU
        %% single GPU
        [outputs_blocks.data,outputs_blocks.id]...
            = worker(fun,Noutputs, use_shared_memory,param,block_id,1,inputs_block); 
    else
        %% run asynchronous on multiple GPUs 
        outputs_blocks(thread_id) = parfeval(@worker,2,fun,Noutputs, use_shared_memory,param,block_id,thread_id,inputs_block); 
    end
end

function out = get_block(full_array, ind, ROI, use_shared_memory)
    Nlayers = length(ind);
    positions = zeros(Nlayers,2,'int32');
    if ~isempty(ROI)
        Npix_small = [length(ROI{1}),length(ROI{2}),length(ind)];
        assert(all(Npix_small <= size(full_array)), 'Provided ROI is larger than input array')
        positions = positions+int32([ROI{1}(1),ROI{2}(1)]-1);
    else
        Npix_small = [size(full_array,1),size(full_array,2),length(ind)];
    end
    % allocate in RAM array to load the block from "stack_object"
    stack_object_block_in = zeros(Npix_small, 'like', full_array);

    if use_shared_memory 
        s = shm(); % create sharemem object 
        if isa(full_array, 'single')
            s.allocate(stack_object_block_in); % allocate SHM memory 
            % attach the shared memory
            [s, stack_object_block_shm] = s.attach(); 
            % ===  write data ===== 
            % use self-made MEX OMP function to move the data 
            utils.get_from_3D_projection(stack_object_block_shm, full_array, positions, ind);
        elseif isa(full_array, 'uint16')           
            stack_object_block_in = zeros(Npix_small, 'like', full_array);
            % ===  write data ===== 
            % use self-made MEX OMP function to move the data 
            utils.get_from_3D_projection(stack_object_block_in, full_array, positions, ind);
            % convert from fp16 to singles
            stack_object_block_in = fp16.get(stack_object_block_in); 

            % upload to the shared memory 
            s.upload(stack_object_block_in); % seems to be really slow ????

            %s.allocate(stack_object_block_in); % allocate SHM memory 
            %[s, stack_object_block_shm] = s.attach(); 
            % % move to shared memory , a bit faster than the memcpy in the
            % % sharedmem function 
            %utils.get_from_3D_projection(stack_object_block_shm, stack_object_block_in, zeros(Nlayers,2,'int32'), int32(1:Nlayers));
            %toc
        end
        % detach the shared memory
        s.detach; 
        % return structure with shared mem object  
        out = s; 
    else
        % use custom made MEX OMP function to move the data 
        stack_object_block_in = utils.get_from_3D_projection(stack_object_block_in, full_array, positions, ind);
        % if needed, convert from fp16 to singles
        out = fp16.get(stack_object_block_in); 
    end
end

function myCleanupFun()
    global status 
    if status
        warning('!! cleaning all shared memory !!')
        % destroy all shared memory that could have been left behind 
        !ipcs -m | cut -d' ' -f2 | grep '^[0-9]' | while read x; do ipcrm -m $x; done
    end
end

function [thread_id, outputs] = gather_block(outputs,param,Nlayers,Noutputs, ind_all, Nl_per_block, output_package)
   if isa(output_package, 'parallel.FevalFuture')
       % gather results from cluster , WAIT FOR CALCULATIONS TO BE FINISHED
%        [~, outputs_block, id] = fetchNext(output_package);
      %% my version of the fetchNext function, it seems faster 
        id = [];
       while true
           for thread_id =1:length(output_package)
               if  strcmpi(output_package(thread_id).State, 'finished') && output_package(thread_id).Read == 0
                  try
                   [outputs_block, id] = output_package(thread_id).fetchOutputs; 
                  catch err
                      if strcmpi(err.identifier, 'parallel:fevalqueue:InvalidExecutionResult')
                          warning('Unknown error, trying to restart parpool')
                          delete(gcp('nocreate'));
                      end
                      if ~isempty(output_package(thread_id).Diary)
                          fprintf('============ THREAD %i FAILED, OUTPUT: ============= \n', thread_id)
                          disp(output_package(thread_id).Diary)
                      end
                      fprintf('============ THREAD %i FAILED, ERROR: ============= \n', thread_id)
                      disp(getReport(output_package(thread_id).Error))
                                 
                      keyboard
                                  
                      output_package.cancel

                      rethrow(err)
                  end
                   break
               end
           end
           if ~isempty(id); break; end
           pause(0.01) % wait for the data to be prepared
       end
   else
      thread_id = 1; 
      id = output_package.id; 
      outputs_block = output_package.data; 
   end
      
   ind = ind_all{id}; 
   
   
   for jj = 1:Noutputs
        % get the results from GPU to RAM 
        % allocate RAM for the final output arrays , allocate when
        % with the same complexity as the outputs_block{ii}{jj} has       
        s = shm();
        if isa(outputs_block{jj}, 'shm')
            % first download from the shared mem back to workspace 
            [s, outputs_block{jj}] = outputs_block{jj}.attach; 
            s.protected = false; 
        end

        % processing returns the same number of layers as input
        is_block = ndims(outputs_block{jj}) == 3 && size(outputs_block{jj},3) == length(ind) ;
        % processing returns one layer per input 
        is_reduced = ismatrix(outputs_block{jj}) &&  size(outputs_block{jj},3) == 1  && size(outputs_block{jj},1) ~= length(ind) &&  ~isempty(param.reduce_fun) && jj ==1; 
        if isa(outputs_block{jj}, 'double')
            % small arrays / scalars just convert
            outputs_block{jj} = single(outputs_block{jj}); 
            % for large complain 
            if is_block
                warning('Output of the function should be blocks in single precision')
            end
        end
        if is_block && param.use_fp16 && isa(outputs_block{jj},'single')
            % move to fp16 precision to avoid large memory allocation 
            outputs_block{jj} = fp16.set(outputs_block{jj}); 
        end

        
        % if the first iteration when it is empty create the large output arrays 
        if isempty(outputs{jj})
            try                           
                if is_block
                    % return large array of size stack_object
                    Npix_small = [size(outputs_block{jj},1),size(outputs_block{jj},2),Nlayers];
                    outputs{jj} = zeros(Npix_small, 'like', outputs_block{jj});
                elseif is_reduced
                    outputs{jj} = [];
                elseif size(outputs_block{jj},1) == length(ind) && ndims(outputs_block{jj}) <= 3
                    % return smaller array split along first dimension 
                    outputs{jj} = zeros([Nlayers, size(outputs_block{jj},2), size(outputs_block{jj},3)], 'like', outputs_block{jj});
                else
                    error('Unimplemented option')
                end
            catch err
               utils.check_available_memory
               rethrow(err)
            end
        end
        
        if is_block
            % if the output is 3D array, add the results into the large output stored in RAM 
            outputs{jj} = tomo.set_to_array(outputs{jj}, outputs_block{jj}, Nl_per_block*(id-1), false);
        elseif is_reduced
            if isempty(outputs{jj})
                outputs{jj} = outputs_block{jj}; 
            else
                outputs{jj} = param.reduce_fun(outputs{jj}, outputs_block{jj}); 
            end
        else
            outputs{jj}(ind,:,:) = outputs_block{jj}; 
        end
   end
end

function [outputs,id] = worker(fun, Noutputs, use_sharemem, param, id, thread_id, inputs)
    
    % let parfor to choose which GPU use
    gpu_id = param.GPU_list(thread_id); 
    gpu = gpuDevice();
    if ~isempty(param.GPU_list) && gpu.Index ~= gpu_id && gpu_id > 0
        gpuDevice(gpu_id);  % avoid unneeded initalization 
    end
        
    if ~isempty(getCurrentTask())
        % report only if inside parfor 
        utils.check_available_memory
        if gpu_id > 0
            fprintf('GPU %i  =====  Memory %3.2GB/%3.2gGB \n',gpu.Index, gpu.AvailableMemory/1e9, gpu.TotalMemory/1e9 )
        end
    end
    
    %% upload data on GPU 

    for jj = 1:length(inputs)
        if isa(inputs{jj}, 'shm')
            % data are downloaded from shared memory 
            [s,inputs{jj}]=inputs{jj}.attach;
        end

        if isnumeric(inputs{jj}) && numel(inputs{jj}) > 1e4 && param.use_GPU && param.move_to_GPU
            inputs{jj} = gpuArray(inputs{jj}); 
        end
    end
    outputs = cell(Noutputs,1); 
    Nlayers = size(inputs{1},3);
    try
        %%%%%% CALL THE WORKING FUNCTION WITH ARRAYS ALREADY MOVED TO GPU %%%% 
        [outputs{:}] = fun(inputs{:});  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if ~isempty(param.reduce_fun)
            % apply additional function to reduce the data 
            switch func2str(param.reduce_fun)
                case 'plus',outputs{1} = sum(outputs{1},3); 
                case 'min', outputs{1} = min(outputs{1},[],3); 
                case 'max', outputs{1} = max(outputs{1},[],3);    
            end
        end
    catch err

        if any(strcmpi(err.identifier, {'parallel:gpu:array:OOMForOperation', 'parallel:gpu:array:OOM'}))
           [~,out] = system(sprintf('nvidia-smi --id=%i | head -n10 | tail -n7', gpu.Index-1)); 
           % disp(getReport(err))
           
           disp(out)
           warning('Low memory on GPU %i, RESETTING GPU ...  \n', gpu.Index)
           
           for jj = 1:length(inputs)
              inputs{jj} = gather(inputs{jj}); 
           end
           % try to clean some leftover data from previous calculations
           reset(gpuDevice)
           
           fprintf('TRYING TO RECURSIVELLY CALL BLOCKFUN %s WITH 2 SUBBLOCKS \n', func2str(fun))
           [outputs{:}] = tomo.block_fun(fun,inputs{:},struct('Nblocks',2));  
        else
            rethrow(err)
        end
    end
    
    
    for jj = 1:length(outputs)
        %% return data from GPU 
        if isa(outputs{jj}, 'gpuArray') && param.use_GPU
            outputs{jj} = gather(outputs{jj}); 
        end
        %% move data to shared memory if needed 
        is_block = ndims(outputs{jj}) == 3 && (size(outputs{jj},3) == Nlayers); 
        if is_block && use_sharemem
            % data are distributed to shared memory 
            s = shm(true);
            s.upload(outputs{jj})
            outputs{jj} = s; 
        end
    end
end
    
function outputs = process_in_single_block(fun, inputs,apply_ROI, splitable_along_3rd_axis, Noutputs, param)
    for jj = 1:length(inputs)
        if splitable_along_3rd_axis(jj)  
            % if stored as fp16 precision, convert to singles 
            inputs{jj} = fp16.get( inputs{jj}); 
            % if provided and possible, use param.ROI to crop the working / inputs array 
            if apply_ROI(jj)
                inputs{jj} = inputs{jj}(param.ROI{:},:);
            end
        end
        if isa(inputs{jj}, 'double')
            inputs{jj} = single(inputs{jj}); 
        end

        if param.use_GPU && isnumeric(inputs{jj}) && numel(inputs{jj}) > 1e2 && param.move_to_GPU
            % upload on GPU 
            inputs{jj} = gpuArray(inputs{jj});
        end
    end
    outputs = cell(Noutputs,1); 
    % execute the worker 
    [outputs{:}] = fun(inputs{:}); 
    if ~isempty(param.reduce_fun)
        % apply additional function to reduce the data 
        switch func2str(param.reduce_fun)
            case 'plus',outputs{1} = sum(outputs{1},3); 
            case 'min', outputs{1} = min(outputs{1},[],3); 
            case 'max', outputs{1} = max(outputs{1},[],3);    
        end
    end
    for jj = 1:Noutputs
        if param.use_GPU && isa(outputs{jj}, 'gpuArray')
            % get back from GPU 
            outputs{jj} = gather(outputs{jj});
        end
        % processing returns the same number of layers as input
        is_block = ndims(outputs{jj}) == 3 && size(outputs{jj},3) == size(inputs{1},3);
        if param.use_fp16 && is_block
            outputs{jj} = fp16.set(outputs{jj}); 
        end
    end
end
