%LOAD_ASTRA_RECONS Load tomographic reconstructions from astra toolbox
% Currently supports reconstruction produced on Summit@OLCF with MPI parallelization.
% Created by YJ.
%

function [recon_astra, par_astra]  = load_astra_recons( par_astra, par, np, Nblocks)

par_astra.dir = strcat('rec_',par_astra.algorithm,'_Niter',num2str(par_astra.Niter));
par_astra.file_name = strcat('rec_',par_astra.algorithm,'_Niter',num2str(par_astra.Niter));

if isfield(par_astra,'min_con') 
    par_astra.dir = strcat(par_astra.dir,'_minCon',num2str(par_astra.min_con));
    par_astra.file_name = strcat(par_astra.file_name,'_minCon',num2str(par_astra.min_con));
end

if isfield(par_astra,'upsample_method') && ~isempty(par_astra.upsample_method)
    par_astra.dir = strcat(par_astra.dir,'_',par_astra.upsample_method);
    par_astra.file_name = strcat(par_astra.file_name,'_',par_astra.upsample_method);
end

par_astra.dir = strcat(par_astra.dir,'_v',num2str(par_astra.N_y),'_h',num2str(par_astra.N_x));


par_astra.time_recon = zeros(Nblocks,1);
par_astra.time_save = zeros(Nblocks,1);

if par_astra.sp>1
    par_astra.dir  = strcat(par_astra.dir,'_sp',num2str(sp));  
    par_astra.file_name = strcat(par_astra.file_name,'_sp',num2str(sp)); 
end
par_astra.file_name = strcat(par_astra.file_name,'_np',num2str(np)); 

recon_temp = {};%

Nslice = 0;
disp('Loading astra reconstructions...')
for i = 1:Nblocks
    disp(strcat(par_astra.file_name,'_',num2str(i),'.hdf5'))
    try
        recon_temp{i} = h5read(strcat(par.output_path,'/tomograms/',par_astra.dir,'/',par_astra.file_name,'_',num2str(i),'.hdf5'),'/rec');
    catch
        disp(strcat(par.output_path,'/tomograms/',par_astra.dir,'/',par_astra.file_name,'_',num2str(i),'.hdf5'))
        
    end
    par_astra.time_recon(i) = sum(h5read(strcat(par.output_path,'/tomograms/',par_astra.dir,'/',par_astra.file_name,'_',num2str(i),'.hdf5'),'/t_recon'));
    par_astra.time_save(i) = sum(h5read(strcat(par.output_path,'/tomograms/',par_astra.dir,'/',par_astra.file_name,'_',num2str(i),'.hdf5'),'/t_save'));
    disp(strcat('recon time:',num2str(par_astra.time_recon(i)/60),'mins.','save time:',num2str(par_astra.time_save(i)/60),'mins'))
    Nslice = Nslice + size(recon_temp{i},3); 
end

disp('Processing astra reconstructions...')

recon_astra = zeros(size(recon_temp{i},2),size(recon_temp{1},1),Nslice,'single');
ind_start = 1;
for i = 1:Nblocks
    disp(strcat(par_astra.file_name,'_',num2str(i),'.hdf5'))
    %ind_start = (i-1)*size(rec,1)+1;
    ind_end = ind_start + size(recon_temp{i},3)-1;
    if strcmp(par_astra.algorithm,'SIRT') || strcmp(par_astra.algorithm,'SART') || strcmp(par_astra.algorithm,'FBP')%2D algorithms
        recon_astra(:,:,ind_start:ind_end) = single(flipud(permute(recon_temp{i},[2 1 3])));
    else
        recon_astra(:,:,ind_start:ind_end) = single((permute(recon_temp{i},[2 1 3])));
    end
    ind_start = ind_end+1;
end

end