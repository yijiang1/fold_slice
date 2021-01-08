%  [U,S,V,rec_all] = SART_SVD(sinogram, theta, Npix, blocks,  par)
%     perform temporal SVD analysis and SART based reconstruction to
%     estimate changes of the sample during reconstruction 
% Inputs: 
%   **sinogram        unwrapped sinogram 
%   **theta          tomography angles 
%   **Npix           size of the reconstructed volume 
%   **blocks         cell list containing indices for each subtomogram
% Outputs: 
%   ++U,S,V          singular vectors 
%   ++rec_all        SVD filterd reconstruction for each subtomogram
% Example of use: 
%   subtomo_ind = [1, find(abs(diff(theta))> 170), length(theta)]; 
%   for ii = 1:length(subtomo_ind)-1
%        ind{ii} = subtomo_ind(ii):subtomo_ind(ii+1); 
%   end
%   [U,S,V] = nonrigid.SART_SVD(sinogram, theta, Npix, ind); 


%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2018 by Paul Scherrer Institute (http://www.psi.ch)    |
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

function [U,S,V,rec_all] = SART_SVD(sinogram, theta, Npix, blocks, varargin)


p = inputParser;
p.addOptional('split', 1)
p.addParameter('valid_angles', [])
p.addParameter('SART_grouping',  25 )  % size of blocks in SART, ART=1, SIRT=Nangles
p.addParameter('GPU', [])   % list of GPUs to be used in reconstruction
p.addParameter('verbose', 1)   % verbose = 0 : quiet, verbose : standard info , verbose = 2: debug 
p.addParameter('N_SVD_modes', 2) % number of recovered SVD modes, 2 is usually enough 
p.addParameter('Niter_SVD', 3)   % number of iter of the SVD SART
p.addParameter('Niter_SART', 5)   % number of internal iterations in each SART loops
p.addParameter('output_folder', '')   % path where the results should be stored 
p.addParameter('mask', [])   % mask applied on the reconstruction 



p.parse(varargin{:})
res = p.Results;


utils.verbose(1,'Using FBP for initial guess')

Nblocks = length(blocks); 

tomogram = cell(Nblocks,1); 
for ii = 1:Nblocks
    utils.progressbar(ii,Nblocks)
       
    % choose projections to process 
    rec_ind = setdiff(blocks{ii}, res.valid_angles); 

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%  

    [Nlayers,width_sinogram,~] = size(sinogram);
    [cfg, vectors] = astra.ASTRA_initialize([Npix,Npix, Nlayers],[Nlayers,width_sinogram],theta); 
    % find optimal split of the dataset for given GPU 
    split = astra.ASTRA_find_optimal_split(cfg, length(res.GPU), 1);

    % new FBP code 
    subtomogram = tomo.FBP_zsplit(sinogram, cfg, vectors, split,'valid_angles',rec_ind,...
        'determine_weights', true, ...
        'GPU', res.GPU ,'filter','ram-lak', 'filter_value',1, 'verbose',-1);

    num_proj_all(ii) = length(rec_ind);

    % get full reconstruction (for FBP is sum already final tomogram)
    % calculate complex refractive index 
    tomogram{ii} = gather(subtomogram);
    
end


if isempty(res.mask)
    constraint_fnct= @(x)x;
else
    constraint_fnct = @(x)(abs(x).*res.mask);
end


for ii = 1:Nblocks
    tomogram{ii} = constraint_fnct(tomogram{ii}); 
end

gpu = gpuDevice; 

for iter = 1:res.Niter_SVD
    utils.verbose(1,' ====== Iteration %i/%i ==== ', iter,res.Niter_SVD)
    rec_all = cat(4, tomogram{:});
    
    utils.verbose(2,'Available GPU memory = %3.1fGB', gpu.AvailableMemory/1e9)
       
    rec_all = reshape(rec_all, [], Nblocks); 
    
    %% %%%%%%%%%%%%%%%%%% APPLY SVD CONSTRAINT %%%%%
    utils.verbose(0,'Calculating SVD ... ')
    Nmodes = min(iter,  res.N_SVD_modes); 
    [U,S,V] = math.fsvd(rec_all,  Nmodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if Nmodes == res.N_SVD_modes
        err_total(iter,:) = gather(sqrt(sum((U*S*V'-rec_all).^2)));
        
        %% plot convergence progress  
        plotting.smart_figure(3)
        loglog(mean(err_total'))
        axis tight
        grid on 
        title('SVD SART - Convergence evolution')
        xlabel('Iteration')
        ylabel('Residuum between SVD model and reconstruction')
        drawnow 

    end


    utils.verbose(0,'Calculating SART ... ')

    % apply SART refinement 
    for ii = 1:Nblocks
        utils.progressbar(ii,Nblocks)


        % choose projections to process 
        rec_ind = setdiff(blocks{ii}, res.valid_angles); 

        if isempty(rec_ind); continue; end

        [cache_SART,cfg_SART] = tomo.SART_prepare(cfg, vectors(rec_ind,:), res.SART_grouping, 'keep_on_GPU', true, 'verbose', 0);

        rec = U*S*V(ii,:)'; 
        rec = reshape(rec,size(tomogram{1})); 

        % get full reconstruction (for FBP is sum already final tomogram)
        % calculate complex refractive index 

        rec = utils.Garray(rec); 
        sino = utils.Garray(sinogram(:,:,rec_ind));
        clear err
        
        for jj = 1:res.Niter_SART
            [rec,err(jj,:)] = tomo.SART(rec, sino, cfg_SART, ...
                vectors(rec_ind,:),cache_SART, 'relax', 0, 'constraint', constraint_fnct, 'verbose', 0);
        end
        % apply some weak total variation to help againts undersampling
        % artefacts 
        %rec = regularization.local_TV3D_chambolle(rec, 1e-6, 10); 
        
        tomogram{ii} = gather(rec); 
        
    end
    clear rec sino cache_SART



end


%% plot SVD evolution 

rec_all = reshape(U*S*V', Npix, Npix,size(sinogram,1), Nblocks); 
rec_all = reshape(rec_all, [size(tomogram{1}), Nblocks]); 

V_sign = sign(mean(V)); 

U(:,1) = U(:,1).*V_sign(1);
V(:,1) = V(:,1).*V_sign(1);

screensize = get( 0, 'Screensize' );


plotting.smart_figure(11)
subplot(1,2,1)
plotting.imagesc3D(squeeze(rec_all(:,:,ceil(end/2),:)))
axis image off 
colormap bone
caxis(gather(math.sp_quantile(rec_all, [0.001, 0.995], 10)));
plotting.suptitle('Tomogram evolution in each subtomogram (central slice)')
subplot(1,2,2)
plot(V, '-o')
title('Principal components evolution')
axis tight 
grid on 
xlabel('Block')
ylabel('S*V''')
Energy = diag(S);
Energy = Energy / sum(Energy); 
for kk = 1:res.N_SVD_modes
   legend_txt{kk} = sprintf('E=%3.2g%%', Energy(kk)*100); 
end
legend(legend_txt ,'location','best')
set(gcf,'Outerposition',[1 screensize(4)-500 800 500]);

if ~isempty(res.output_folder) && ~debug()
   try
       savefig(fullfile(res.output_folder, 'SVD_filtered_evolution.fig'))
   catch err 
      warning('Saving of SVD_filtered_evolution failed with error: %s', err.message)
   end
end

U = reshape(U, [size(tomogram{1}), res.N_SVD_modes]); 

U_plot = U(:,:,2:end-1,:);  % it seems that first and last layer are not well estimated  
U_plot = U_plot - median(quantile(min(U_plot,[],1),0.01,2),3); 
U_plot = U_plot ./ median(quantile(max(U_plot,[],1),0.99,2),3); 

% 
U_plot = cat(2,  U_plot(:,:,:,1), U_plot(:,:,:,2)); 



plotting.smart_figure(10)
subplot(2,1,1)
plotting.imagesc3D(U_plot, 'init_frame', size(U_plot,3)/2)
axis image off 
colormap bone
caxis(gather(math.sp_quantile(U_plot, [0.001, 0.995], 10)));
title('Principal components (left is 1th PC , right is 2nd PC)')
subplot(2,1,2)
plotting.imagesc3D(U_plot, 'init_frame', size(U_plot,1)/2,  'slider_axis',1)
axis image off 
colormap bone
title('Principal components (left is 1th PC , right is 2nd PC)')
caxis(gather(math.sp_quantile(U_plot, [0.001, 0.995], 10)));

%%%suptitle('Principal vectors showing tomogram evolution (slide to see layers of the sample)')
set(gcf,'Outerposition',[1 screensize(4)-1250 1200 700]);
if ~isempty(res.output_folder)  && ~debug()
    print('-f10','-dpng','-r300',[res.output_folder, '/SVD_modes_scaled.png']);
end



%% get reconstructions to RAM 
U = gather(U);
S = gather(S); 
V = gather(V); 
rec_all = gather(rec_all); 

U = reshape(U, [Npix, Npix,size(sinogram,1),res.N_SVD_modes]); 


end