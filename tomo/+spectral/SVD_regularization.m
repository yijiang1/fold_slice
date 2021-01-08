%  [U,S,V,rec_all, rec_all_0] = SVD_regularization(sinogram, theta, Niter_SVD, Npix, reconstruct_ind,max_projections,par)
%     perform spectral SVD analysis and SART based reconstruction 
% Inputs: 
%   sinogram - unwrapped sinogram 
%   theta - projection angles 
%   Niter_SVD - number of iteration of the SVD optimization 
%   Npix - (3x1 int), size of the output volume 
%   reconstruct_ind - (int array) list of angles that will be considered for reconstruction 
%   max_projections - maximal number of projections selected for each
%       energy step. Set "inf" to ignore this limit. It is useful to
%       equialize the weight for each energy step when the distribution is
%       highly unequal 
%   par - tomogrpahy paramter structure 
% Outputs: 
%  U,S,V - SVD vectors 
%  rec_all - regularized SART reconstructions for each energy 
%  rec_all_0 - oriignal FBP reocnstructions before SVD regularization


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


function [U,S,V,rec_all, rec_all_0] = SVD_regularization(sinogram, theta, Niter_SVD, Npix, reconstruct_ind,max_projections,par)

% internal parameters 
SART_grouping = 25; 
Niter_SART = 10; 
Nmodes = 2; 
Nangles = length(theta); 

E_all = unique(par.energy(ismember(1:Nangles, reconstruct_ind) )); 
Nenergy = length(E_all);
tomogram_edensity = cell(Nenergy,1); 

for ii = 1:Nenergy
    utils.progressbar(ii,Nenergy)
       
    % choose projections to process 
    rec_ind = find(par.energy == E_all(ii) & ismember(1:Nangles,reconstruct_ind)');   % use only some angles 
    % downsample to the requested "max_projections"
    rec_ind = rec_ind(1:max(1,ceil(length(rec_ind)/max_projections)):end); 
    
    num_proj_all(ii) = length(rec_ind); 
    if isempty(rec_ind)
        tomogram_edensity{ii} = zeros(Npix,Npix,size(sinogram,1),'single'); 
        continue; 
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%  

    [Nlayers,width_sinogram,~] = size(sinogram);
    CoR = [Nlayers,width_sinogram]/2;   % there is 0.5px shift between fft_1d/fft_2d vs none
    [cfg, vectors] = astra.ASTRA_initialize([Npix,Npix, Nlayers],[Nlayers,width_sinogram],theta,par.lamino_angle,par.tilt_angle,1,CoR); 
    % find optimal split of the dataset for given GPU 
    split = astra.ASTRA_find_optimal_split(cfg, length(par.GPU_list), 1);

    % new FBP code 
    tomogram = tomo.FBP_zsplit(sinogram, cfg, vectors, split,'valid_angles',rec_ind,...
        'determine_weights', true, ...
        'GPU', par.GPU_list,'filter',par.filter_type, 'filter_value',par.freq_scale, 'verbose',0);


    % Caclulate delta tomogram 

    par.lambda = 1.234e-9 / E_all(ii);

    par.factor=par.lambda/(2*pi*par.pixel_size); 
    par.factor_edensity = 1e-30*2*pi/(par.lambda^2*2.81794e-15);

    % get full reconstruction (for FBP is sum already final tomogram)
    % calculate complex refractive index 
    tomogram_edensity{ii} = gather(tomogram*par.factor*par.factor_edensity);
    
end

% quantity = min(max_projections, hist(par.energy, E_all));

rec_all_0 = gather(cat(4, tomogram_edensity{:}));

% get a weighted average 
mrec = max(0,mean(rec_all_0 .* reshape(num_proj_all,1,1,1,[]) ,4)) / mean(num_proj_all); 
mask = mrec > graythresh(mrec);

%% 
mask = imopen(mask, strel('disk',5)); 
mask = imdilate(mask, strel('disk',5)); 
mask = imfill(mask, 'holes');
mask = Garray(mask);
constraint_fnct = @(x)(max(0,x.*mask));

for ii = 1:Nenergy
    tomogram_edensity{ii} = constraint_fnct(tomogram_edensity{ii}); 
end

gpu = gpuDevice; 

for iter = 1:Niter_SVD
    utils.verbose(1,' ====== Iteration %i ==== ', iter)
    rec_all = cat(4, tomogram_edensity{:});
    
    utils.verbose(2,'Available GPU memory = %3.1fGB', gpu.AvailableMemory/1e9)
    
        
    plotting.smart_figure(4)
    plotting.imagesc3D(cat(2, squeeze(rec_all(:,:,ceil(end/2),:)), squeeze(rec_all_0(:,:,ceil(end/2),:)))); 
    axis off image, colormap bone
    caxis(gather(math.sp_quantile(rec_all, [0.001, 0.9999], 10)));
    title('Left - SVD refined , Right - original FBP')

    
    
    rec_all = reshape(rec_all, [], Nenergy); 
    
    % add more modes progressivelly for higher iteration number 
    Nmodes_tmp = min(floor(iter^(1/3)), Nmodes); 
    
    
    %%%%%%%%%%%%%%%%%%%% APPLY SVD CONSTRAINT %%%%%
    [U,S,V] = math.fsvd(rec_all, Nmodes_tmp);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% enforce smoothness of the first V vector 
   
    
    if iter > 1
        % plot convergence progress  
        try
            err_total(iter,:) = gather(sqrt(sum((U*S*V'-rec_all).^2)));
        catch
           keyboard 
        end
        plotting.smart_figure(3)
        loglog(mean(err_total'))
        axis tight
        grid on 
        xlabel('Iteration')
        ylabel('Residuum between SVD model and reconstruction')
        drawnow 
    end

    
    
   
    
%     relax = 0.1; 
%     V(:,1) = V(:,1) *(1-relax)  + relax * polyval(polyfit(1:Nenergy,V(:,1)',1), 1:Nenergy)';
%  
%     if Nmodes_tmp == 3
%         keyboard
%        
%     end



    if iter > 1
        % plot SVD evolution 
            
        best_angle = fminsearch(@(x)get_rotation_score(x,U,S,V), zeros(3,1)); 
        R = rotation_matrix_3D(best_angle(1),best_angle(2),best_angle(3)); 
        R = R(1:Nmodes_tmp,1:Nmodes_tmp); % use only number of modes that is needed
        Vplot = (R*V')'; % apply rotation 
        Uplot = U*S*inv(R);  % apply inverse rotation  

        % flip sign to for convinince 
        sign_U = sign(mean(Uplot)); 
        Uplot = Uplot .* sign_U;
        Vplot = Vplot .* sign_U;
        rec_all = gather(rec_all);

        U_plot = Uplot ./ quantile(Uplot(1:32:end,:), 0.99); 
        
        
        U_plot = reshape(U_plot, [size(tomogram_edensity{1}),Nmodes_tmp]); 
        U_plot = reshape(permute(U_plot,[1,2,4,3]), [size(U_plot,1),size(U_plot,2)*size(U_plot,4),size(U_plot,3)]); 
       
        plotting.smart_figure(5)
        subplot(2,1,1)
  
        plotting.imagesc3D(U_plot, 'init_frame', size(U_plot,3)/2)
        axis image off 
        colormap bone
        caxis(gather(math.sp_quantile(U_plot, [0.001, 0.995], 10)));
        title('Topos SVD vectors')
        subplot(2,1,2)
        plot(E_all,Vplot, '-o')
        hold all        
        offset = min(min(Vplot )); 
        range = max(max(Vplot )) - offset; 
        bar(E_all,num_proj_all / max(num_proj_all)*range*0.2+ offset, 'facecolor', 'none', 'Basevalue', offset)
        hold off 
        title('Chronos SVD vectors')
        axis tight 
        grid on 
        xlabel('Energy [keV]')
        ylabel('S*V''')
        Energy = diag(S);
        Energy = Energy / sum(Energy); 
        for kk = 1:Nmodes_tmp
           legend_txt{kk} = sprintf('E=%3.2g%%', Energy(kk)*100); 
        end
        legend(legend_txt ,'location','best')
        clear U_plot
    end


    % ather from GPU to save memory 
%     Uout = gather(Uout);    
    

    % apply SART refinement 
    for ii = 1:Nenergy
        utils.progressbar(ii,Nenergy)


        % choose projections to process 
        rec_ind = find(par.energy == E_all(ii) & ismember(1:Nangles,reconstruct_ind)');   % use only some angles 
        % downsample to the requested "max_projections"
        rec_ind = rec_ind(1:ceil(length(rec_ind)/max_projections):end); 

        if isempty(rec_ind); continue; end
        

        [cache_SART,cfg_SART] = tomo.SART_prepare(cfg, vectors(rec_ind,:), SART_grouping, 'keep_on_GPU', true);

        
        rec = U*S*V(ii,:)'; 
        rec = reshape(rec,size(tomogram_edensity{1})); 

        par.lambda = 1.234e-9 / E_all(ii);
        par.factor=par.lambda/(2*pi*par.pixel_size); 
        par.factor_edensity = 1e-30*2*pi/(par.lambda^2*2.81794e-15);

        % get full reconstruction (for FBP is sum already final tomogram)
        % calculate complex refractive index 
        rec = rec / (par.factor*par.factor_edensity);

        rec = Garray(rec); 
        sino = Garray(sinogram(:,:,rec_ind));
        clear err
        
        for jj = 1:Niter_SART
            [rec,err(jj,:)] = tomo.SART(rec, sino, cfg_SART, ...
                vectors(rec_ind,:),cache_SART, 'relax', 0, 'constraint', constraint_fnct);
        end
        % apply some weak total variation to help agsints undersampling
        % artefacts 
        rec = regularization.local_TV3D_chambolle(rec, 1e-6, 10); 

        tomogram_edensity{ii} = (rec*par.factor*par.factor_edensity);

%         plotting.smart_figure(1)
%         subplot(1,3,1)
%         plotting.imagesc3D(tomogram_edensity{ii}); axis off image, colormap bone 
%         caxis([0,1])
%         subplot(1,3,2)
%         plot(theta(rec_ind), 'o-')
%         title(['Nangles ', num2str(length(rec_ind)) ])
%         subplot(1,3,3)
%         loglog(mean(err'))
%         drawnow 
    end
    clear rec sino cache_SART
end

rec_all = reshape(rec_all, [size(tomogram_edensity{1}), Nenergy]); 




% get recosntructions to RAM 
U = gather(U);
S = gather(S); 
V = gather(V); 
rec_all = gather(rec_all); 

end

function score = get_rotation_score(x,U,S,V)
    R = rotation_matrix_3D(x(1),x(2),x(3)); 
    Nmodes = size(S,1); 
    V = (R(1:Nmodes,1:Nmodes)*V')';
    score = gather(norm(V(:,1) - mean(V(:,1)))); 
end
