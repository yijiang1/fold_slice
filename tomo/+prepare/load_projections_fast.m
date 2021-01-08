%   LOAD_PROJECTIONS_FAST load reconstructed projections from disk to RAM 
%
%  [stack_object, theta,num_proj, par]  =  load_projections_fast(par, exclude_scans, dims_ob, theta, custom_preprocess_fun)
%
%   Inputs: 
%       **par - parameter structure 
%       **exclude_scans - list of scans to be excluded from loading, [] = none
%       **dims_ob  - dimension of the object 
%       **theta - angles of the scans 
%       **custom_preprocess_fun - function to be applied on the loaded data, eg cropping , rotation, etc 
%
%   *returns* 
%       ++stack_object - loaded complex-valued projections 
%       ++theta - angles corresponding to the loaded projections, angles for missing projections are removed 
%       ++num_proj - number of projections 
%       ++par - updated parameter structure 


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
%   in written form in the publication: "Data processing was carried out 
%   using the "cSAXS matlab package" developed by the CXS group,
%   Paul Scherrer Institut, Switzerland." 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided "as they are" without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.



function [stack_object, theta,num_proj, par]  = load_projections_fast(par, exclude_scans, dims_ob, theta, custom_preprocess_fun)

import ptycho.* 
import utils.* 
import io.*
import plotting.imagesc3D
utils.verbose(struct('prefix', 'loading'))

if nargin < 5
    custom_preprocess_fun = []; 
end
if ~isempty(custom_preprocess_fun)  && ishandle(custom_preprocess_fun) && ~strcmpi(func2str(custom_preprocess_fun), '@(x)x')
    custom_preprocess_fun = [] ;
end

scanstomo = par.scanstomo; 

% avoid loading scans listed in 'exclude_scans'
if  ~isempty(exclude_scans)
    ind = ismember(scanstomo, exclude_scans); 
    %% clear values corresponding to excluded scans 
    scanstomo(ind) = []; 
    theta(ind) = []; 
    par.subtomos(ind) = [];
end
    


% % plot average vibrations for each of the loaded projections 
% utils.verbose(0,'Checking stability of the projections')
% poor_projections = prepare.plot_sample_stability(par, scanstomo, ~par.online_tomo, par.pixel_size); 
% if sum(poor_projections) && ...
%         (par.online_tomo || ~strcmpi(input(sprintf('Remove %i low stability projections: [Y/n]\n',sum(poor_projections)), 's'), 'n') )
%     theta(poor_projections) = []; 
%     scanstomo(poor_projections) = []; 
% else
%     utils.verbose(0,'All projections are fine')
% end



verbose(0,'Checking available files')
verbose(0); % make it quiet 
missing_scans = []; 
proj_file_names = cell(length(scanstomo),1);
for num = 1:length(scanstomo)
    progressbar(num, length(scanstomo))
    filename = find_projection_files_names(par, scanstomo(num));
    if isempty(filename)
        missing_scans(end+1) = scanstomo(num); 
        continue
    end
    proj_file_names{num} = filename;

end

verbose(par.verbose_level); % return to original settings

if ~isempty(missing_scans)
    plotting.smart_figure(1)
    subplot(2,1,1)
    hold on 
    plot(missing_scans, theta(ismember(scanstomo, missing_scans)), 'rx', 'Linewidth', 2)
    hold off
    legend({'Measured angles', 'Missing scans'})
    axis tight
    drawnow 
    
    ind = ismember(scanstomo, missing_scans); 
    verbose(1,['Scans  not found are ' num2str(missing_scans)])
    verbose(1,['Projections not found are ' num2str(find(ind))])
    %% clear values corresponding to measured but missing scans (not reconstructed) 
    scanstomo(ind) = []; 
    theta(ind) = []; 
    proj_file_names(ind) = []; 
    par.subtomos(ind) = [];
else
    verbose(1,'All projections found')
end

num_proj = length(scanstomo); 

if isfield(par, 'fp16_precision') && par.fp16_precision
    dtype = uint16(1i);  % use uint32 to store half floar precision data
else
    dtype= single(1i);
end

downsample = 2^par.downsample_projections; % calculate downsample factor for binning , default par.downsample_projections = 0; 

pixel_size =zeros(num_proj,2);
energy = zeros(num_proj,1); 
stack_object=zeros(ceil(dims_ob(1) / downsample),ceil(dims_ob(2) / downsample),num_proj, 'like', dtype);
residua = zeros(num_proj,1); 

%disp(size(stack_object))
tic


% load at least 10 frames per worker to use well the resources 
block_size =  max(1, feature('numcores'))*4; 

object_ROI = {ceil(1+par.asize(1)/2/downsample):ceil((dims_ob(1)-par.asize(1)/2)/downsample),ceil(1+par.asize(2)/2/downsample):ceil((dims_ob(2)-par.asize(2)/2)/downsample)}; 

verbose(1,'Loading projections ...')

missing_all = []; 
t0 = tic;
%% load data, use parfor but process blockwise to avoid large memory use and also allow user stopping during MEX reading 
for block_id = 1:ceil(num_proj/block_size)
    block_inds = 1+(block_id-1)*block_size: min(num_proj, block_id*block_size); 
    utils.progressbar(block_id, ceil(num_proj/block_size))
   
     if strcmpi(par.file_extension, 'h5')  && ~verLessThan('matlab', '9.4') && ... 
                    (isfield(par, 'use_mex_loader') && par.use_mex_loader ) % only matlab newer than R2018a is supported
         object_block = []; 
         try
             % fast MEX loader, sometimes it tends to fail and needs to
             % be run again to load the data corectly
            [object_block,missing_tmp] = mex_read(par.dims_ob_loaded,  proj_file_names(block_inds), par.Nthreads_mexread); 
         catch Err
            disp('Error in loading using MEX, falling back to matlab reader, try to reduce par.Nthreads_mexread is this warning repeats often')
            disp(Err)
         end

         % if loading was not succeful .. 
         if isempty(object_block)
            [object_block, missing_tmp] = matlab_read(par.dims_ob_loaded,  proj_file_names(block_inds)); 
         end
    else
        % loading using matlab for original MAT file data or old matlab 
        [object_block, missing_tmp] = matlab_read(par.dims_ob_loaded,  proj_file_names(block_inds)); 
    end
    missing_all = [missing_all, block_inds(missing_tmp)];

        
    % read additional information  
    for jj = setdiff(block_inds, block_inds(missing_tmp)) % remove missing projection from loading
        if strcmpi(par.file_extension, 'h5')
            try
            pixel_size(jj,:) = h5read(proj_file_names{jj}, '/reconstruction/p/dx_spec');
            energy(jj) = h5read(proj_file_names{jj}, '/reconstruction/p/energy'); 
            catch err
                disp(err)
                keyboard
            end 
        else
            % load it from the matlab file is not supported (it is too slow)
            pixel_size(jj,:) = par.pixel_size; 
            energy(jj) = nan; 
        end
    end
        
    %% apply custom data proprocessing and caculate basic statistics, DO IT ON GPU 
    [object_block, residua(block_inds,1), projection_value(block_inds)] = ...
        tomo.block_fun(@process_projection_block, object_block, custom_preprocess_fun, par, object_ROI,pixel_size(block_inds,:), struct('verbose_level', 0)); 

    % convert data to fp16 precision if requested 
    if isfield(par, 'fp16_precision') && par.fp16_precision
        object_block = fp16.set(object_block); 
    end
    
    
    
    stack_object(:,:,block_inds) = object_block; 
    
end

pixel_size = min(pixel_size,[],2); % projection were already rescaled to provide the same pixel size in each dimension 

verbose(1, 'Data loaded in %is', round(toc(t0)))

% downsample the illum_sum if requested 
if downsample > 0 
    par.illum_sum = utils.binning_2D(crop_pad(par.illum_sum, ceil(dims_ob/downsample)*downsample) , downsample); 
    par.asize = ceil(par.asize / downsample); 
end




failed_projections = projection_value < 0.1*median(projection_value) | ismember(1:num_proj, missing_all) | ~isfinite(projection_value); 

if any(failed_projections )
    verbose(0,['Projections failed are ' num2str(find(failed_projections))])
    verbose(0,['Scans failed are       ' num2str(scanstomo(failed_projections))])
else
    verbose(0,'All loaded projections seems OK')
end


[Nprojections] = size(stack_object,3); 
poor_projections = false;

if ~par.is_laminography
    % laminography has a more complex definition of field of view ->
    % currently not implemented 
    
    poor_projections = (residua(:)' > par.max_residua_limit) ;   % ignore in the case of laminography 
    verbose(1, 'Found %i/%i projections with more than %i residues ', sum(poor_projections), Nprojections, par.max_residua_limit)
    if any(poor_projections)
        verbose(1,['Projections with residua are ' num2str(find(poor_projections))])
        verbose(1,['Scans with residua are       ' num2str(scanstomo(poor_projections))])
    end
end
verbose(1, 'Find residua done')


% avoid also empty projections
which_remove =  poor_projections | failed_projections; 

%%% Getting rid of missing projections %%%
if any(which_remove)
    
    [Nx,Ny,~] = size(stack_object);
    title_extra = {}; 
    for ii = 1:num_proj
        if which_remove(ii)
            title_extra{end+1} = sprintf(' N residua: %i',residua(ii)); 
        end
    end
    
    verbose(1,' %i failed projections shown in figure(1) \n', sum(which_remove))

    tomo.show_projections(stack_object(:,:,which_remove), theta(which_remove), par, 'fnct', @angle, ...
        'title', 'Projection to be removed','plot_residua', true, 'title_extra', title_extra,  ...
        'rectangle_pos', [par.asize(2)/2,Ny-par.asize(2)/2, par.asize(1)/2,Nx-par.asize(1)/2], 'figure_id', 1) 


    
    if par.online_tomo || debug() || ~strcmpi(input(sprintf('Do you want remove %i failed/wrong projections and keep going (Y/n)?',sum(which_remove)),'s'),'n') 
        verbose(0,'Removing failed/wrong projections. stack_object, scanstomo, theta and num_proj are modified')
        
        stack_object(:,:,which_remove) = [];
        scanstomo(which_remove)=[];
        theta(which_remove)=[];
        pixel_size(which_remove,:) = []; 
        energy(which_remove,:) = []; 
        residua(which_remove)=[];
        par.subtomos(which_remove) = [];
        
        verbose(0,'Done')
    else
        verbose(0,'Keeping empty spaces for failed projections. Problems are expected if you continue.')
    end
end


par.scanstomo = scanstomo; 
par.num_proj=numel(scanstomo);

% store number of residua for later processing
par.nresidua_per_frame = residua(:)'; 

pixel_scale = pixel_size ./ min(pixel_size(:)); % just in case that the axis do not have indetical pixel size, NOT TESTED YET 

assert(par.num_proj > 0, 'No projections loaded')


if all(all(abs(pixel_scale)-1 < 1e-6)) || ~any(isfinite(mean(pixel_scale)))
    %if all datasets have the same pixel scale 
    pixel_scale = [1,1]; 
else
    warning('Datasets do not have equal pixel sizes, auto-rescaling projections')
    % use FFT base rescaling -> apply illumination function first to remove
    % effect of the noise out of the reconstruction region
    rot_fun = @(x,sx,sy)(utils.imrescale_frft(x .*  par.illum_sum, sx, sy)) ./ ( max(0,utils.imrescale_frft(par.illum_sum,sx,sy))+1e-2*max(par.illum_sum(:))); 
    stack_object = tomo.block_fun(rot_fun,stack_object, pixel_scale(:,1),pixel_scale(:,2));
    pixel_scale = [1,1]; 
end


par.pixel_scale = pixel_scale; 
par.energy = energy; 



if size(stack_object,3) ~= length(theta) || length(theta) ~= par.num_proj
    error('Inconsistency between number of angles and projections')
end
utils.verbose(struct('prefix', 'template'))


end

function [object_block, residua, projection_value] = process_projection_block(object_block, custom_preprocess_fun, par, object_ROI, pixel_size)
    % auxiliary function used to apply various preprocessing steps, ie custom_preprocess_fun, binning, clipping and residua calculation on the
    % object_block on GPU -> avoid CPU-GPU transfer overhead 
    % returns: 
    %   object_block - processed complex valued projections 
    %   residua - number of residua in each frame 
    %   projection_value - average amplitude of the projection 
    %   pixel_size in each dimension 

    % apply additional processing, e.g. rotation 
    if ~isempty(custom_preprocess_fun) 
        object_block = custom_preprocess_fun(object_block); 
    end
    
    if any(pixel_size(:,1) ~= pixel_size(:,2))
        % in the case of asymmetric pixel size, 
        % upsample the data in the dimennsion with lower resolution (-> at least relax issues in tomography interpolation)
        pixel_scale = pixel_size ./ min(pixel_size,[],2) ;
        assert(all(std(pixel_scale) < 1e-3), 'Variable resolution between projection and asymmetric pixel size is not implemented')
    
        Npix = size(object_block); 
        dims_ob_new = round(Npix(1:2) .* pixel_scale(1,:)); 
        object_block = utils.interpolateFT(object_block, dims_ob_new);
    end
    
    Npix = size(object_block); 
    downsample = 2^par.downsample_projections; 
    % downsample the data if requested 
    if downsample > 1
        object_block = utils.binning_2D(utils.crop_pad(object_block, ceil(Npix/downsample)*downsample) , downsample); 
    end
    
    %% clip the projections amplitude by quantile filter 
    if par.clip_amplitude_quantile > 0 && par.clip_amplitude_quantile < 1
        MAX = quantile(reshape(abs(object_block(1:10:end,1:10:end,:)), [], Npix(3)), par.clip_amplitude_quantile ,1);
        MAX = reshape(MAX,1,1,[]);
        clip_fun = @(x,M)(min(abs(x),M) .* x ./ (abs(x) + 1e-5));
        object_block = clip_fun(object_block, MAX);
    end

    residua = squeeze(math.sum2(abs(utils.findresidues(object_block(object_ROI{:},:)))>0.1));
    projection_value = squeeze(math.sum2(abs(object_block)));

end

function [object_block, missing] = matlab_read(dims_ob,  proj_file_names)
    % projection loading using matlab 
    % Inputs: 
    %  dims_ob - projection size 
    %  proj_file_names - cell of filenames  to be loaded 
    % Outputs: 
    %  object_block - loaded projection
    %  missing - list of missing (failed) projections 
    
    object_block = zeros([dims_ob,length(proj_file_names)], 'like', single(1i));  
    loaded = false(length(proj_file_names),1); 
    for jj = 1:length(proj_file_names)
        utils.verbose(2,['Reading file: ' proj_file_names{jj}])
        try
            object = io.load_ptycho_recons(proj_file_names{jj}, 'object');
            object = single(object.object); 
            object = prod(object,4);   % use only the eDOF object if multiple layers are available
            object_block(:,:,jj) = utils.crop_pad(object, dims_ob); 
            loaded(jj) = true;
        catch
            utils.verbose(-1,'Loading of file %s failed', proj_file_names{jj})
        end
    end
    missing = find(~loaded); 
end



function [object_block, missing] = mex_read(dims_ob,  proj_file_names, Nthreads)
    % fast projection loader by MEX with paralelization 
    % Inputs: 
    %  dims_ob - projection size 
    %  proj_file_names - cell of filenames  to be loaded 
    %  Nthreads - number of threads used to load the projections in parallel
    % Outputs: 
    %  object_block - loaded projection
    %  missing - list of missing (failed) projections 
    
    Nthreads = min(length(proj_file_names),Nthreads ); 
    
    proj_file_names  = reshape(proj_file_names, 1,[]); 
    
    % load complex-valued projections using parallel MEX 
    try
        [object_block, missing] = io.ptycho_read(Nthreads, 'single', dims_ob, '/reconstruction/object', proj_file_names); 
    catch err 
        if strcmpi(err.identifier, 'ptycho:read:failed')
            Nthreads = 5; 
            [object_block, missing] = io.ptycho_read(Nthreads, 'single', dims_ob, '/reconstruction/object', proj_file_names); 
            warning off backtrace
            warning('===================================================================================================================================')
            warning('Loading of projections failed due to too high multithreading, if this warning repeats, consider lowering par.Nthreads_mexread value')
            warning('===================================================================================================================================')
            warning on backtrace
        else
           rethrow(err) 
        end
    end
    assert(ndims(object_block) <= 5, 'Unexpected dimensionality of inputs')
    
    object_block = permute(object_block, [2,1,3,4,5]); % transpose loaded reconstructions 
    object_block = prod(object_block,4) ; % get one eDoF frame if ML reconstruction is used 
    object_block = squeeze(object_block); % get rid of extra dimensions 

    assert(ndims(object_block) == 3, 'Unexpected dimensionality of inputs')
    
 end

