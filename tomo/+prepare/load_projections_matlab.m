%   LOAD_PROJECTIONS_MATLAB load reconstructed projections from disk to RAM 
%  created by YJ based on PSI's function
%  [stack_object, theta,num_proj, par]  =  load_projections_matlab(par, exclude_scans, dims_ob, theta, custom_preprocess_fun)
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

function [stack_object, theta,num_proj, par]  = load_projections_matlab(par, exclude_scans, dims_ob, theta, custom_preprocess_fun)

import ptycho.* 
import utils.* 
import io.*
import plotting.*

if nargin < 5
    custom_preprocess_fun = []; 
end
if ~isempty(custom_preprocess_fun)  && ishandle(custom_preprocess_fun) && ~strcmpi(func2str(custom_preprocess_fun), '@(x)x')
    custom_preprocess_fun = [] ;
end

scanstomo = par.scanstomo; 
if isfield(par,'energy')
    energy = par.energy; 
else
    energy = zeros(length(theta),1); 
end

% avoid loading scans listed in 'exclude_scans'
if  ~isempty(exclude_scans)
    ind = ismember(scanstomo, exclude_scans); 
    scanstomo(ind) = [];
    theta(ind) = [];
    energy(ind) = [];
end
    
% % plot average vibrations for each of the laoded projections 
% disp('Checking stability of the projections')
% poor_projections = prepare.plot_sample_stability(par, scanstomo, ~par.online_tomo, par.pixel_size); 
% if sum(poor_projections) && ...
%         (par.online_tomo || ~strcmpi(input(sprintf('Remove %i low stability projections: [Y/n]\n',sum(poor_projections)), 's'), 'n') )
%     theta(poor_projections) = []; 
%     scanstomo(poor_projections) = []; 
% else
%     disp('All projections are fine')
% end

verbose(1,'Checking available files')
missing_scans = []; 
proj_file_names = {};
proj_recon_method = {};
proj_roi = {};
proj_scanNo = {};
for num = 1:length(scanstomo)
    progressbar(num, length(scanstomo))
    %proj_file_names{num} = find_ptycho_filename(par.analysis_path,scanstomo(num),par.fileprefix,par.filesuffix, par.file_extension);
    %proj_file_names{num} = find_projection_files_names_aps(par, scanstomo(num));
    [proj_file_names{num},proj_recon_method{num},proj_roi{num},proj_scanNo{num}] = find_ML_recon_files_names(par, scanstomo(num));
    %disp(proj_file_names{num})
    if isempty(proj_file_names{num})
        missing_scans(end+1) = scanstomo(num); 
    end
end

verbose(par.verbose_level); % return to original settings
%{
figure(1)
subplot(2,1,1)
hold on 
plot(missing_scans, theta(ismember(scanstomo, missing_scans)), 'rx')
hold off
legend({'Measured angles', 'Missing projections'})
axis tight
%}

if ~isempty(missing_scans)
    ind = ismember(scanstomo, missing_scans); 
    verbose(1,['Scans  not found are ' num2str(missing_scans)])
    verbose(1,['Projections not found are ' num2str(find(ind))])
    scanstomo(ind) = []; 
    theta(ind) = []; 
    proj_file_names(ind) = []; 
    proj_recon_method(ind) = [];
    proj_roi(ind) = [];
    proj_scanNo(ind) = [];
    energy(ind) = [];
else
    verbose(1,'All projections found')
end

num_proj = length(scanstomo); 
object_size_orig = zeros(2,num_proj);
if isfield(par, 'fp16_precision') && par.fp16_precision
    % use uint32 to store half floar precision data
    stack_object=zeros(dims_ob(1),dims_ob(2),num_proj, 'like', fp16.set(1i));
else
    stack_object=zeros(dims_ob(1),dims_ob(2),num_proj, 'like', single(1i));
end
pixel_scale =zeros(num_proj,2);

tic

if num_proj == 0
    verbose(0, 'No new projections loaded')
    return
end

which_missing = false(1,num_proj);     % Include here INDEX numbers that you want to exclude (bad reconstructions)

utils.check_available_memory
%%
wb = waitbar(0,'1','Name','Loading ptycho-tomo projection...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(wb,'canceling',0);
%
t0 = tic;

for num=1:num_proj
    % Update waitbar and message
    status = sprintf(par.scan_string_format, scanstomo(num));
    status = strcat(status,' (',num2str(num),'/',num2str(num_proj),') ');
    
    if num>1
        timeLeft = (num_proj-num+1)*avgTimePerIter;

        if timeLeft>3600
            time_status = sprintf(' Time left:%3.3g hour', timeLeft/3600);
        elseif timeLeft>60
            time_status = sprintf(' Time left:%3.3g min', timeLeft/60);
        else
            time_status = sprintf(' Time left:%3.3g sec', timeLeft);
        end
        status = strcat(status,time_status);
    end
    waitbar(num/num_proj,wb,status)
    
    % Check for clicked Cancel button
    if getappdata(wb,'canceling')
        break
    end
    file = proj_file_names{num};

    if ismember(scanstomo(num), exclude_scans)
        warning(['Skipping by user request: ' file{1}])
        continue  % skip the frames that are listed in exclude_scans
    end
    
    if ~iscell(file)
        file = {file};  % make them all cells 
    end
    object= [];
    for jj = length(file):-1:1
        for ii=1:3
            try
                object = load(file{1},'object');
                object = single(object.object);

                parameter = load(file{1},'p');
                pixel_scale(num,:) = parameter.p.dx_spec;  %pixel size
                %disp(file{1})
                break
            catch
                warning(['Loading failed: ' [file{1}]])
            end
        end
    end
    
    %% for multislice recon - sum layers into a single projection
    if size(object,3)>1
        if isfield(par.MLrecon,'select_layers') && any(par.MLrecon.select_layers)
            object = prod(object(:,:,par.MLrecon.select_layers),3);
        else
            object = prod(object,3);
        end
    end
    
    %%
    object_size_orig(:,num) = size(object);

    if isempty(object) || all(object(:) == 0 )
        which_missing(num) = true;
        warning(['Loading failed: ' [file{:}]])
        continue
    end
    if isfield(par, 'crop_edge') && par.crop_edge>0
        object = object(1+par.crop_edge:end-par.crop_edge,1+par.crop_edge:end-par.crop_edge);
    end
    if ~isempty(custom_preprocess_fun) 
        object = custom_preprocess_fun(object); 
    end
    
    nx = dims_ob(2);
    ny = dims_ob(1);

    if size(object,2) > nx       
        object = object(:,1:nx);       
    elseif size(object,2) < nx
        object = padarray(object,[0 nx-size(object,2)],'post');
    end

    if size(object,1) > ny
        if par.auto_alignment|| par.get_auto_calibration
            object = object(1:ny,:);
        else
            shifty = floor((size(object,1)-ny)/2);
            object = object([1:ny]+shifty,:);
        end
    elseif size(object,1) < ny
        if par.auto_alignment||par.get_auto_calibration
            object = padarray(object,[ny-size(object,1) 0],'post');
        else
            shifty = (ny-size(object,1))/2;
            object = padarray(object,[ny-size(object,1)-floor(shifty) 0],'post');
            object = padarray(object,[floor(shifty) 0],'pre');
        end
    end
    
    stack_object(:,:,num) = object;

    % if par.showrecons
    %     mag=a+bs(object);
    %     phase=angle(object);
    %     figure(1); clf
    %     imagesc(mag); axis xy equal tight ; colormap bone(256); colorbar; 
    %     title(['object magnitude S',sprintf('%05d',ii),', Projection ' ,sprintf('%03d',num) , ', Theta = ' sprintf('%.2f',theta(num)), ' degrees']);drawnow;
    %     set(gcf,'Outerposition',[601 424 600 600])
    %     figure(2); imagesc(phase); axis xy equal tight; colormap bone(256); colorbar; 
    %     title(['object phase S',sprintf('%05d',ii),', Projection ' ,sprintf('%03d',num) , ', Theta = ' sprintf('%.2f',theta(num)), ' degrees']);drawnow;
    %     set(gcf,'Outerposition',[1 424 600 600])    %[left, bottom, width, height
    %     figure(3); %  imagesc3D(probe); 
    %     axis xy equal tight
    %     set(gcf,'Outerposition',[600 49 375 375])    %[left, bottom, width, height
    %     figure(4); 
    %     if isfield(p, 'err')
    %         loglog(p.err);
    %     elseif isfield(p, 'mlerror') 
    %         loglog(p.mlerror)
    %     elseif isfield(p, 'error_metric')
    %         loglog(p.error_metric(2).iteration,p.error_metric(2).value)
    %     end
    %     title(sprintf('Error %03d',num))
    %     set(gcf,'Outerposition',[1 49 600 375])    %[left, bottom, width, height
    %     drawnow;
    % end

 	avgTimePerIter = toc(t0)/num;

end  % enf of parfor
delete(wb)
%store info for ML reconstructions
par.proj_file_names = proj_file_names;
par.proj_recon_method = proj_recon_method;
par.proj_roi = proj_roi;
par.proj_scanNo = proj_scanNo;
par.object_size_orig = object_size_orig;
verbose(1, 'Data loaded')

%% examine projections
verbose(1, 'Find residua')
[Nx, Ny, Nprojections] = size(stack_object); 

object_ROI = {ceil(1+par.asize(1)/2:Nx-par.asize(1)/2),ceil(1+par.asize(2)/2:Ny-par.asize(2)/2)}; 
residua = tomo.block_fun(@(x)(squeeze(math.sum2(abs(utils.findresidues(x))>0.1))),stack_object, struct('ROI', {object_ROI})); 

if isfield(par,'max_residua_limit')
   max_residua = par.max_residua_limit; 
else
    max_residua = 100;
end
poor_projections = (residua(:)' > max_residua) & ~par.is_laminography ;   % ignore in the case of laminography 

if any(poor_projections) 
    verbose(1, 'Found %i/%i projections with more than %i residues ', sum(poor_projections), Nprojections, max_residua)
end


if any(which_missing & ~ismember(scanstomo,  exclude_scans) )
    missing = find(which_missing & ~ismember(scanstomo,  exclude_scans)); 
    verbose(1,['Projections not found are ' num2str(missing)])
    verbose(1,['Scans not found are       ' num2str(scanstomo(missing))])
else
    verbose(1,'All projections loaded')
end
toc

% avoid also empty projections
which_wrong =  poor_projections | squeeze(math.sum2(stack_object)==0)'; 

if any(which_wrong & ~ismember(scanstomo,  exclude_scans) )
    wrong = find(which_wrong & ~ismember(scanstomo,  exclude_scans)); 
    verbose(1,['Projections failed are ' num2str(wrong)])
    verbose(1,['Scans failed are       ' num2str(scanstomo(wrong))])
else
    verbose(1,'All loaded projections are OK')
end


%%% Getting rid of missing projections %%%
which_remove = which_missing | which_wrong; 
if any(which_remove)
    if par.online_tomo || ~strcmpi(input(sprintf('Do you want remove %i missing/wrong projections and keep going (Y/n)?',sum(which_remove)),'s'),'n') 
        disp('Removing missing/wrong projections. stack_object, scanstomo, theta and num_proj are modified')
        
        stack_object(:,:,which_remove) = [];
        scanstomo(which_remove)=[];
        theta(which_remove)=[];
        pixel_scale(which_remove,:) = []; 
        energy(which_remove,:) = []; 

        disp('Done')
    else
        disp('Keeping empty spaces for missing projections. Problems are expected if you continue.')
    end
end

par.scanstomo = scanstomo; 
par.num_proj=numel(scanstomo);

pixel_scale = pixel_scale ./ mean(pixel_scale); 

assert(par.num_proj > 0, 'No projections loaded')

if all(all(abs(pixel_scale)-1 < 1e-6)) || ~any(isfinite(mean(pixel_scale)))
    %if all datasets have the same pixel scale 
    pixel_scale = [1,1]; 
else
    warning('Datasets do not have equal pixel sizes!')
    %warning('Datasets do not have equal pixel sizes, auto-rescaling projections')
    
    % use FFT base rescaling -> apply illumination function first to remove
    % effect of the noise out of the reconstruction region
    
    %rot_fun = @(x,sx,sy)(utils.imrescale_frft(x .*  par.illum_sum, sx, sy)) ./ ( max(0,utils.imrescale_frft(par.illum_sum,sx,sy))+1e-2*max(par.illum_sum(:))); 
    %stack_object = tomo.block_fun(rot_fun,stack_object, pixel_scale(:,1),pixel_scale(:,2));
    %pixel_scale = [1,1]; 
end

par.pixel_scale = pixel_scale; 
par.energy = energy; 

%% clip the projections ampltitude by quantile filter 
if par.clip_amplitude_quantile < 1
    MAX = quantile(reshape(abs(fp16.get(stack_object(1:10:end,1:10:end,:))), [], par.num_proj), par.clip_amplitude_quantile ,1); 
    MAX = reshape(MAX,1,1,par.num_proj); 
    clip_fun = @(x,M)(min(abs(x),M) .* x ./ (abs(x) + 1e-5)); 
    stack_object = tomo.block_fun(clip_fun,stack_object, MAX, struct('use_GPU', true));
end

if size(stack_object,3) ~= length(theta) || length(theta) ~= par.num_proj
    error('Inconsistency between number of angles and projections')
end

%{
if ~isempty(par.tomo_id) && all(par.tomo_id > 0)
    % sanity safety check, all loaded angles correpont to the stored angles
    [~,theta_test] = prepare.load_angles(par, par.scanstomo, [], false);
    if max(abs(theta - theta_test)) > 180/par.num_proj/2
        error('Some angles have angles different from expected')
    end
end
%}

%% replot angle
plot_angles = true;
if par.verbose_level && plot_angles
    plotting.smart_figure(1);
    subplot(2,1,1)
    plot(par.scanstomo,theta,'ob'); grid on; 
    xlim(par.scanstomo([1,end]))
    %legend('Tilt angles')
    xlabel('Scan #')
    ylabel('Tilt angles')
    
    %[anglessort,~] = sort(theta);

    subplot(2,1,2)
    plot(diff(theta))
    ylabel('Angle increment')
    %title('Angular spacing'); 
    grid on; 
    xlim([1,par.num_proj-1])
    if par.windowautopos
        screensize = get( groot, 'Screensize' );
        win_size = [946 815]; 
        set(gcf,'Outerposition',[139 min(163,screensize(4)-win_size(2))  win_size]);  %[left, bottom, width, height]
    end
    title('Measured angles')
    drawnow 
end


end
