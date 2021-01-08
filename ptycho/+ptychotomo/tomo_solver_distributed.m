% TOMO_SOLVER_DISTRIBUTED tomography solver that loads reconstruction from
% ptychography, process them and use them to update the current tomgoraphy
% reconstruction. Then create new estimates of ptychography objects 
%
% volData = tomo_solver_distributed(par, volData, projData_rec, p0, theta_all, scan_ids)
% 
%  Inputs: 
%   **par          - parameter structure 
%   **volData      - initial guess of the tomographic volumes, values are linearized, ie projection = exp(sum(volData,1)) = prod(exp(volData),1)
%   **projData_rec - initial guess of the projections for each angle, probe, positions, etc 
%   **p0           - basic parameters for the ptychography solver
%   **theta_all    - all angles in degrees 
%   **scan_ids     - scan numbers that correspons to each angle value 
%  *returns*
%   ++volData      - final reconstruction of the tomographic volume 

% Academic License Agreement
%
% Source Code
%
% Introduction 
% •	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
%   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the PtychoShelves 
%   computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").
%
% Terms and Conditions of the LICENSE
% 1.	LICENSOR grants to LICENSEE a royalty-free, non-exclusive license to use the PROGRAM for academic, non-commercial purposes, upon the terms and conditions 
%       hereinafter set out and until termination of this license as set forth below.
% 2.	LICENSEE acknowledges that the PROGRAM is a research tool still in the development stage. The PROGRAM is provided without any related services, improvements 
%       or warranties from LICENSOR and that the LICENSE is entered into in order to enable others to utilize the PROGRAM in their academic activities. It is the 
%       LICENSEE’s responsibility to ensure its proper use and the correctness of the results.”
% 3.	THE PROGRAM IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR 
%       A PARTICULAR PURPOSE AND NONINFRINGEMENT OF ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. IN NO EVENT SHALL THE LICENSOR, THE AUTHORS OR THE COPYRIGHT 
%       HOLDERS BE LIABLE FOR ANY CLAIM, DIRECT, INDIRECT OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION WITH THE PROGRAM OR THE USE 
%       OF THE PROGRAM OR OTHER DEALINGS IN THE PROGRAM.
% 4.	LICENSEE agrees that it will use the PROGRAM and any modifications, improvements, or derivatives of PROGRAM that LICENSEE may create (collectively, 
%       "IMPROVEMENTS") solely for academic, non-commercial purposes and that any copy of PROGRAM or derivatives thereof shall be distributed only under the same 
%       license as PROGRAM. The terms "academic, non-commercial", as used in this Agreement, mean academic or other scholarly research which (a) is not undertaken for 
%       profit, or (b) is not intended to produce works, services, or data for commercial use, or (c) is neither conducted, nor funded, by a person or an entity engaged 
%       in the commercial use, application or exploitation of works similar to the PROGRAM.
% 5.	LICENSEE agrees that it shall make the following acknowledgement in any publication resulting from the use of the PROGRAM or any translation of the code into 
%       another computing language:
%       "Data processing was carried out using the PtychoShelves package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
%       Scherrer Institut, Switzerland."
%
% Additionally, any publication using the package, or any translation of the code into another computing language should cite 
% K. Wakonig, H.-C. Stadler, M. Odstrčil, E.H.R. Tsai, A. Diaz, M. Holler, I. Usov, J. Raabe, A. Menzel, M. Guizar-Sicairos, PtychoShelves, a versatile 
% high-level framework for high-performance analysis of ptychographic data, J. Appl. Cryst. 53(2) (2020). (doi: 10.1107/S1600576720001776)
% and for difference map:
% P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379–382 (2008). 
%   (doi: 10.1126/science.1158573),
% for maximum likelihood:
% P. Thibault and M. Guizar-Sicairos, Maximum-likelihood refinement for coherent diffractive imaging, New J. Phys. 14, 063004 (2012). 
%   (doi: 10.1088/1367-2630/14/6/063004),
% for LSQ-ML:
% M. Odstrčil, A. Menzel, and M. Guizar-Sicairos, Iterative least-squares solver for generalized maximum-likelihood ptychography, Opt. Express 26(3), 3108 (2018). 
%   (doi: 10.1364/OE.26.003108),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
%   (doi: 10.1364/OE.24.029089),
% and/or for OPRP:
% M. Odstrcil, P. Baksh, S. A. Boden, R. Card, J. E. Chad, J. G. Frey, W. S. Brocklesby,  Ptychographic coherent diffractive imaging with orthogonal probe relaxation. 
% Opt. Express 24.8 (8360-8369) 2016. (doi: 10.1364/OE.24.008360).
% 6.	Except for the above-mentioned acknowledgment, LICENSEE shall not use the PROGRAM title or the names or logos of LICENSOR, nor any adaptation thereof, nor the 
%       names of any of its employees or laboratories, in any advertising, promotional or sales material without prior written consent obtained from LICENSOR in each case.
% 7.	Ownership of all rights, including copyright in the PROGRAM and in any material associated therewith, shall at all times remain with LICENSOR, and LICENSEE 
%       agrees to preserve same. LICENSEE agrees not to use any portion of the PROGRAM or of any IMPROVEMENTS in any machine-readable form outside the PROGRAM, nor to 
%       make any copies except for its internal use, without prior written consent of LICENSOR. LICENSEE agrees to place the following copyright notice on any such copies: 
%       © All rights reserved. PAUL SCHERRER INSTITUT, Switzerland, Laboratory for Macromolecules and Bioimaging, 2017. 
% 8.	The LICENSE shall not be construed to confer any rights upon LICENSEE by implication or otherwise except as specifically set forth herein.
% 9.	DISCLAIMER: LICENSEE shall be aware that Phase Focus Limited of Sheffield, UK has an international portfolio of patents and pending applications which relate 
%       to ptychography and that the PROGRAM may be capable of being used in circumstances which may fall within the claims of one or more of the Phase Focus patents, 
%       in particular of patent with international application number PCT/GB2005/001464. The LICENSOR explicitly declares not to indemnify the users of the software 
%       in case Phase Focus or any other third party will open a legal action against the LICENSEE due to the use of the program.
% 10.	This Agreement shall be governed by the material laws of Switzerland and any dispute arising out of this Agreement or use of the PROGRAM shall be brought before 
%       the courts of Zürich, Switzerland.


function volData = tomo_solver_distributed(par, volData, projData_rec, p0, theta_all, scan_ids)
    import utils.* 
    
    % be sure to clean GPU first 
    reset(gpuDevice)
    
    p0.queue.path = par.queue_path; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%  create the queue %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    system(['rm -rf ', par.queue_path]);
    % if exist(par.queue_path, 'dir')
    %     rmdir(par.queue_path, 's'); 
    % end
    Nangles = length(theta_all); 

    try; mkdir([par.queue_path, '/pending/']); end
    try; mkdir([par.queue_path, '/done/']); end

    system('rm temp/*/prepared_data.mat ');
    

    % minimum p structure needed to create the queue
    p_init.prepare_data_path = par.prepare_data_path;
    
    % create file queue
    for ii = 1:Nangles
        utils.progressbar(ii,length(scan_ids)); 
        p = p_init; 
        p.scan_number = scan_ids(ii); 
        save('-v6',fullfile(par.queue_path, sprintf('pending/scan%05d.mat',scan_ids(ii))), 'p'); 
    end
    
    
    initial_settings = 'p_initial.mat'; 

    save(initial_settings, 'p0');

    
    
    Npx_vol = size(volData);
    assert(Npx_vol(1)==Npx_vol(2), 'Horizontal volume size has to be symmetric')
    assert(mod(Npx_vol(1),64) == 0, 'Input volum esize  should be dividable by 64')
    thickness = p0.thickness;                                     % Total thickness of the simulated sample if multiple layers are provided; 

    Nangles = length(projData_rec); 
    
    % choosed optimally distributed the angles in the 180 deg range to minimize overlap 
    angle_step = median(diff(theta_all)); 
    init_angle = math.argmin(abs(theta_all - 180-angle_step*par.downsample_angles/2));
    angle_indices =   [1:par.downsample_angles:init_angle-1, init_angle:par.downsample_angles:Nangles];
    volData_hist = {[],[]};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    
    % initialize one dataset
    [p0] = core.initialize_ptycho(p0);
    p0.prepare_data_path = par.prepare_data_path;

    
    gpu = gpuDevice();

    
    
    %% create some loost support mask r
    if par.apply_support
        m_volData = abs(mean(volData,3)); 
        support_mask = utils.imgaussfilt2_fft(m_volData, 20) > mean(m_volData(:))/2; 
        support_mask = utils.imgaussfilt3_conv(support_mask, [10,10,0]);
        volData = volData .* support_mask; 
    else
        support_mask = 1;  
    end
    %% move arrays to GPU 
    par.support_mask = Garray(support_mask); 
    par.norm_full = norm(volData(:)); 

    % move to uint8 - the volume has to be on GPU !! moving from and to GPU
    % causes too much overhead 
    volData = ptychotomo.compress_volume(volData, par);    
    volData = Garray(volData); 
    
    % prepate support masks for the projections using the "support_mask" calculated from the
    % reconstructed tomography volume 
    par.weight = sum(support_mask); 
    par.weight = utils.crop_pad(par.weight, [1,p0.object_size(2)]); 
    par.weight = 1-par.weight / max(par.weight); 
    par.support_mask = support_mask; 
    
    % store all volumes per iteration 
    volData_all = cell(par.Niter,1); 
    
            

    disp('Reconstructing')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% iterativelly solve ptycho tomo task %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    update_norm     = nan(par.Niter, Nangles, par.Niter_inner);
    fourier_error   = nan(par.Niter, Nangles,2);
    update_norm_acc = nan(par.Niter, 1);
    hist_numbers    = zeros(Nangles,1);
    ptycho_results = cell(Nangles,1); % store temporally data used for reconstruction , avoid data storing  
    

    Npx_proj = [p0.object_size, 1];

    proj_number_list = []; % list of projections yet in procesing

    
    t_plot = tic; 
    for iter = 1:par.Niter
        if iter >= max(par.ptycho_reconstruct_start, par.ptycho_ML_reconstruct_start)
            Nlayers = min(2^ceil((iter+1 - par.ptycho_ML_reconstruct_start)), par.Nlayers_max ); 
        else
            Nlayers = 1; 
        end

        Npx_proj(3) = Nlayers;


        full_reconstruction = mod(iter,par.ptycho_interval) == 0 && iter >= par.ptycho_reconstruct_start; 
          
        fprintf('Free GPU mem: %3.4gGB \n', gpu.AvailableMemory / 1e9)
        
        fprintf('Iteration %i/%i Nlayers %i full_reconstruction %i \n', iter, par.Niter, Nlayers, full_reconstruction)
        layer_distance = thickness / Nlayers * ones(1,Nlayers-1) ; 
        
        Ngroups = length(angle_indices); 

        group_ind = get_golden_ratio_groups(theta_all(angle_indices), Ngroups);

                
        %%%%% start accelerated gradients %%%%%%%%%%%%%%
        if iter == par.ptycho_accel_start
            volData_hist{1} = ptychotomo.decompress_volume(gather(volData), par); 
            volData_hist{2} = volData_hist{1}; 
        elseif iter > par.ptycho_accel_start
            beta = (iter-par.ptycho_accel_start+1)/(iter-par.ptycho_accel_start+3);

            volData_hist{1} = volData_hist{2}; 
            volData_hist{2} = ptychotomo.decompress_volume(gather(volData), par); 
            
            upd = beta*(volData_hist{1}-  volData_hist{2}); 
            update_norm_acc(iter) = norm(upd(:));

            if ~isfinite(update_norm_acc(iter))
                keyboard
            end
            if all(isfinite(update_norm_acc([iter, iter-1]))) && (update_norm_acc(iter) > update_norm_acc(iter-1))
               % reset the accelerated gradients 
               warning('Reset accelerated gradients')
                volData_hist{1} = ptychotomo.decompress_volume(gather(volData), par); 
                volData_hist{2} = volData_hist{1}; 
                par.ptycho_accel_start = iter; 
            else
                % else move in the accelerated direction 
                volData = volData_hist{2} + upd;
            end
            
            volData = ptychotomo.compress_volume(volData, par); 
            volData = utils.Garray(volData);  % move it back to GPU after gathering for acceleratin computation 
            clear upd
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        verbose(struct('prefix', {'ptychotomo'}))

        for outer_loop = 1:par.Niter_inner
            fprintf('Outer loop %i \n', outer_loop)
            
            tic
        
          group_id = 1; 
          while group_id <= Ngroups || ~isempty(proj_number_list)
            if verbose() < 0 
                utils.progressbar(group_id, Ngroups)
            end
            clear volData_upd_full

            
            verbose(0,'PREPARING')
            try
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%% prepare new initial guess if needed %%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            proj_number_list_update = []; 
            while length(proj_number_list) < par.max_queue_length && group_id <= Ngroups
                         
                ind = angle_indices(group_ind{group_id}); 

                % calculate how the projection should look like and store the initial guess to "initial_guess.mat"
                for ll = 1:length(ind)
                    proj_number = ind(ll); 
                    scan_number = scan_ids(proj_number);

                    if ~isempty(dir(sprintf('%s/pending/scan%05i.mat', par.queue_path,scan_number)))
                        verbose(0, 'Preparing %i', scan_number)

                        verbose(0); 
                        [projData_model{proj_number}, projData_rec{proj_number}, ptycho_results{proj_number}] = ...
                            ptychotomo.prepare_distributed_data(p0, volData, projData_rec{proj_number},layer_distance, par, ...
                                                                iter == par.ptycho_reconstruct_start,  iter >= par.ptycho_reconstruct_start  );
                        
                        if iter >= par.ptycho_reconstruct_start
                            % file is ready for processing, move it to queue
                            io.movefile_fast(sprintf('%s/pending/scan%05i.mat', par.queue_path,scan_number), sprintf('%s/scan%05i.mat', par.queue_path,scan_number));
                        else
                            % skip ptychography reconstrution, pretent that the file is already finished 
                            io.movefile_fast(sprintf('%s/pending/scan%05i.mat', par.queue_path,scan_number), sprintf('%s/done/scan%05i.mat', par.queue_path,scan_number));
                        end
                        hist_numbers(proj_number) = hist_numbers(proj_number) +1; 
                        proj_number_list_update(end+1) = proj_number; 
                        proj_number_list(end+1) = proj_number;
                        continue
                    else
                        verbose(0, 'Waiting for scan %i to be returned to /pending/', scan_number)
                    end
                end
                group_id = group_id + 1; 

            end
            catch err 
                keyboard
            end
           
         
            clear scan_number  proj_number
            
            %%%%%%% run the reconstruction using the generated queue files %%%
            
           
            if iter < par.ptycho_reconstruct_start  
                verbose(0,'TOMO RECONSTRUCTING')
                for proj_number = proj_number_list_update
                    % save previous object estimation and pretend that they were just calculated 
                    proj = projData_rec{proj_number};
                    pout.object{1} = exp(proj.object_c);
                    pout.probes = proj.probe;
                    pout.positions = proj.positions;
                    pout.asize = p0.asize;
                    ferr = nan;
                    ptycho_results{proj_number} = pout; 
                end
                prepare_data_path = sprintf(p0.prepare_data_path, scan_ids(proj_number)); 
            else
                verbose(0,'PTYCHOTOMO RECONSTRUCTING')
            end


            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            verbose(struct('prefix', {'ptychotomo'}))
            verbose(0,'GATHERING')

            %%%%%%%%%%%%%%   apply the recent reconstructions into the 3D volume %%%%%%
            % load the stored reconstructions and return corresponding object update 
            

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% start gathering reconstructions %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            volData_upd_full = 0; 
            proj_number = proj_number_list(1); 
            scan_number = scan_ids(proj_number); 
            proj_number_list(1) = []; 
            scan_file_path = sprintf('%s/done/scan%05d.mat', par.queue_path, scan_number); 


            timeout_gathering = 0.1; % [s]
            ii = 0; 
            while true 
                verbose(0)
                if ii * timeout_gathering > par.wait_time_solver 
                    % skip projetions that were not delivered in more than
                    % "wait_time_solver" time 
                    warning('Failed to gather %s', scan_file_path)
                    break
                end
                if exist(scan_file_path, 'file')
                    break
                else
                    % wait for solver to prepare the reconstruction 
                    pause(timeout_gathering)
                    if ii == 0;  verbose(0,'Waiting for reconstruction %i', scan_number); end
                    ii = ii + 1;
                    continue
                end

            end
            
            if ii * timeout_gathering > par.wait_time_solver 
                 continue
            end
            
            p0.scan_number = scan_number;
            
            
            verbose(0)
            verbose(0, 'Gathering %i', p0.scan_number)

            if iter >= par.start_3D_reconstruction
                par.update_step = par.lambda/length(ind)/outer_loop; 
                par.update_step = complex(par.update_step / 5, par.update_step); % much slower convergence for amplitude 
            else
                par.update_step = 0; 
            end
            
             try
                if ii * timeout_gathering < par.wait_time_solver 
                    % calculate and accumulate update of the 3D volume 
                    [volData,projData_rec{proj_number}, fourier_error(iter,proj_number,:), update_norm(iter,proj_number,outer_loop)] = ...
                            ptychotomo.gather_distributed_reconstructions(volData, projData_model{proj_number}, ptycho_results{proj_number},par);
                end
            catch err 
                warning('gather_distributed_reconstructions failed: %s', err.message)
            end


 


            % file is already used, move it back to pending ... 
            io.movefile_fast(scan_file_path, sprintf('%s/pending/scan%05i.mat', par.queue_path,p0.scan_number));

            
            % release shared memory 
            if isa(ptycho_results{proj_number}, 'shm')
                ptycho_results{proj_number}.protected = false; % make it possible to delete the SHM 
                ptycho_results{proj_number}.detach; 
                ptycho_results{proj_number}.free()  ; % free shared memory         
                ptycho_results{proj_number} = []; 
            end
            
            % save memory 
            projData_model{proj_number}.object = gather(projData_model{proj_number}.object); 
            projData_model{proj_number}.object_c = gather(projData_model{proj_number}.object_c); 
            ptycho_results{proj_number} = []; 

            
            if math.norm2(projData_rec{proj_number}.probe) > 2.1
                keyboard
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% %%%%%%%%%%% PLOT PROGRESS %%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if par.debug ||  toc(t_plot) >  par.plot_every
                plotting.smart_figure(24)
                subplot(2,2,1)
                plotting.imagesc3D(-imag( ptychotomo.decompress_volume(volData(:,:,ceil(end/2)),par))); 
                caxis(gather(math.sp_quantile(-imag(volData(:,:,100)), [5e-3, 1-5e-3], 10)))
                axis off image
                colormap bone 
                title('Current reconstruction')
                subplot(2,2,2)
                ind_ok = find(any(~isnan(fourier_error(:,:,1)),2)); 
                plot(ind_ok,fourier_error(ind_ok,:,1), '-')
                hold all 
                plot(ind_ok,fourier_error(ind_ok,:,2), '--')
                plot(ind_ok,nanmean(fourier_error(ind_ok,:,1),2), 'k-', 'Linewidth', 2)
                plot(ind_ok,nanmean(fourier_error(ind_ok,:,2),2), 'k--', 'Linewidth', 2)
                hold off 
                set(gca, 'xscale', 'log')
                set(gca, 'yscale', 'log')
                grid on 
                axis tight 
                title('Fourier error (from ptycho)')
                subplot(2,2,3)
                plot(update_norm(:,:,1))
                hold all
                plot(update_norm_acc, 'k--', 'Linewidth', 2)
                plot(nanmean(update_norm(:,:,1),2), 'k', 'Linewidth', 2)
                plot(nanmean(update_norm(:,:,end),2), 'k:', 'Linewidth', 2)
                hold off 
                set(gca, 'xscale', 'log')
                set(gca, 'yscale', 'log')
                grid on 
                title('Update difference (in tomo)')
                axis tight 
                subplot(2,2,4)
                plotting.imagesc3D(exp(crop_pad(sum(projData_rec{proj_number}.object_c,3), [Npx_vol(3), Npx_vol(1)])))
                axis xy off image  
                title('Example of a projection')
                %drawnow
                
                plotting.smart_figure(15222)
                subplot(2,2,1)
                plotting.imagesc3D(abs(exp(sum(projData_rec{proj_number}.object_c,3))))
                colorbar
                colormap bone 
                axis xy off image  
                caxis([0.8, 1.1])
                title('Ptycho - abs')
                subplot(2,2,2)
                plotting.imagesc3D(angle(exp(sum(projData_rec{proj_number}.object_c,3))))
                colorbar
                colormap bone 
                axis xy off image  
                title('Model - phase')

                subplot(2,2,3)
                plotting.imagesc3D(abs(exp(sum(projData_model{proj_number}.object_c,3))))
                colorbar
                colormap bone 
                axis xy off image  
                title('Model - abs')
                caxis([0.8, 1.1])
                subplot(2,2,4)
                plotting.imagesc3D(angle(exp(sum(projData_model{proj_number}.object_c,3))))
                colorbar
                colormap bone 
                title('Model - phase')

                axis xy off image  
                
                
                     
                if  par.Niter_inner > 1
                    plotting.smart_figure(13123)
                    plot(squeeze(update_norm(max(1,iter-1),:,:))')
                    hold all 
                    plot(nanmedian(squeeze(update_norm(max(1,iter-1),:,:)),1) ,'k','Linewidth', 2) 
                    hold off 
                    set(gca, 'xscale', 'log')
                    set(gca, 'yscale', 'log')
                    title(sprintf('Inner loop convergence iter = %i Nlayers = %i', iter, Nlayers))
                    axis tight 
                end
                
                
                drawnow 
                            
                t_plot = tic; 
            end
            

          end
        
        toc
        


        
        end
        
        try


        volData = ptychotomo.decompress_volume(volData, par); 
        
        % plotting.smart_figure(123123)
        % subplot(1,2,1)
        % plotting.imagesc3D(log(mean(abs(math.fftshift_2D(fft2(exp(gather(volData))))),3)))
        % axis image off
        % colormap(plotting.colormaps.franzmap)
        % title('FFT of reconstruction')
        % subplot(1,2,2)
        % plotting.imagesc3D(log(mean(abs(math.fftshift_2D(fft2(exp(gather(volData_upd_full))))),3)))
        % axis image off 
        % colormap(plotting.colormaps.franzmap)
        % title('FFT of last update')
        % suptitle(sprintf('Iter %i Nlayers %i fullrecons %i',  iter, Nlayers, full_reconstruction))
        
        
        %%
        plotting.smart_figure(234)
        frame = ceil(size(volData,3)/2); 
        subplot(1,2,1)
        plotting.imagesc3D(-imag(volData(:,:,frame)));
        caxis(gather(math.sp_quantile(-imag(volData(:,:,frame)), [5e-3, 1-5e-3], 10)))
        axis off image 
        colorbar
        colormap bone 
        title('Phase')
        subplot(1,2,2)
        plotting.imagesc3D(-real(volData(:,:,frame)));
        caxis(gather(math.sp_quantile(-real(volData(:,:,frame)), [5e-3, 1-5e-3], 10)))
        colorbar
        axis off image 
        colormap bone 
        title('Absorbtion')
        plotting.suptitle(sprintf('Reconstruction in %i-th iteration', iter))
        %%
        
        img = -gather(imag(volData(:,:,frame))); 
        img = img / math.sp_quantile(img, 1-5e-3, 1); 
        img = max(0, min(1, img)); 
        if ~exist('./results_ptychotomo', 'dir'); mkdir('results_ptychotomo'); end 
        imwrite(gray2ind(img), bone, sprintf('results_ptychotomo/img_3d_ptycho_iter%i_Nlayers_%i.png', iter, Nlayers))

        
        if ~par.debug
            %% save a preview on disk 
           % suptitle(sprintf('Iter %i Nlayers %i fullrecons %i',  iter, Nlayers, full_reconstruction))
           % print('-dpng', '-f234', '-r300', sprintf('results_ptychotomo/progress_3d_ptycho_iter%i_Nlayers_%i.png', iter, Nlayers))
        end
        
        catch err 
           warning('Plotting failed') 
           keyboard
        end



        volData_all{iter} = gather(volData(:,:,ceil(end/2))); 
                  
        % save memory
        volData_upd_full = []; 

        volData = ptychotomo.compress_volume(volData, par); 
            
    end
        
  

    %% show final results 
    
    for iter = 1:par.Niter-1
        try
        volupd(:,:,iter) =volData_all{iter+1} - volData_all{iter}; 
        volprev(:,:,iter) =volData_all{iter}; 
        end
    end
    
    for iter = 1:par.Niter-2
        try
        a = volupd(:,:,iter);
        b = volupd(:,:,iter+1);
        C(iter)  = abs(corr(a(:), b(:))); 
        end
    end


    figure  
    subplot(1,2,1)
    plotting.imagesc3D(-imag(volprev))
    colormap bone 
    axis off image 
    title('Phase')
    subplot(1,2,2)
    plotting.imagesc3D(-real(volprev))
    colormap bone 
    axis off image 
    title('Absorbtion')
    plotting.suptitle('Reconstruction quality in each iteration')
    
    figure    
    plotting.imagesc3D(volupd)
    axis off image 
    title('Complex volume update in each iteration')

    
    
    
end



