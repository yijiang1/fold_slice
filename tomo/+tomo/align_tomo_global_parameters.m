%    ALIGN_TOMO_GLOBAL_PARAMETERS   find center of rotation or lamino angle or tilt of the projections 
%   plot various statistics that may (and may not) help to decided which
%   parameter provides best reconstruction
%
%   align_tomo_global_parameters(sinogram,angles, Npix, par, varargin )
%
% Inputs: 
%       **sinogram_0        - real value sinogram (ie not diff)
%       **angles            - angle in degress
%       **Npix              - size of the reconstructed field 
%       **par               - parameter structure -> params, INPUTS DESCRIBED IN CODE 
% Outputs: 
%       (none)
%       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!
%       updates should be done manually by user if one is confident that
%       the newly estimated geometry is definitelly leading to improved
%       reconstruction 
%       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%       it very useful for quick verification that the global geometry is ok

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


function  align_tomo_global_parameters(sinogram,angles, Npix, params, varargin )
    


    import tomo.*
    import utils.*
    import math.*
    utils.verbose(struct('prefix', 'align'))


    parser = inputParser;
    parser.addParameter('binning',  4 , @isint )
    parser.addParameter('deformation_fields', [])  % assume deformated sample and use these fielresid_sino
    parser.addParameter('plot_results',  true , @islogical ) % plot results 
    parser.addParameter('verbose',  1 , @isnumeric ) % change verbosity of the code 
    parser.addParameter('is_laminography',  false , @isnumeric ) % change verbosity of the code 
    parser.addParameter('search_range',  [-100,100] , @isnumeric ) % search range for the center of rotation 
    parser.addParameter('num_grid_points',  100 , @isnumeric )    % number of grid points 
    parser.addParameter('search_parameter',  'center_of_rotation' , @(x)(ismember(lower(x), {'center_of_rotation', 'center_of_rotation_y', 'lamino_angle', 'tilt_angle', 'rot_angle', 'shear_angle' })) )
    parser.addParameter('CoR_offset', 0, @isnumeric); 
    parser.addParameter('CoR_offset_v', 0, @isnumeric); 
    parser.addParameter('lamino_angle_offset', 0, @isnumeric); 
    parser.addParameter('tilt_angle_offset', 0, @isnumeric); 
    parser.addParameter('rotation_angle_offset', 0, @isnumeric); 
    parser.addParameter('shear_angle_offset', 0, @isnumeric); 
    parser.addParameter('selected_roi', {}, @iscell); 
    parser.addParameter('usecircle', false, @islogical); 
    parser.addParameter('showed_layer_id', [], @isnumeric); 
    parser.KeepUnmatched = false;
    parser.parse(varargin{:})
    r = parser.Results;
    
    % load all to the param structure 
    par = params; 
    for name = fieldnames(r)'
        if ~isfield(par, name{1}) || ~ismember(name, parser.UsingDefaults) % prefer values in param structure 
            par.(name{1}) = r.(name{1});
        end
    end

    % load all to the param structure 

    verbose(1,'Starting %s estimation', r.search_parameter)
      
    verbose(1,['Binning: ', num2str(r.binning)])

       
    sinogram = tomo.block_fun(@imreduce,sinogram,r.selected_roi,r.binning);
    
    %% %%%%%%%%%%%%%%%% initialize astra %%%%%%%%%%%%%%%%
    [Nlayers,width_sinogram,~] = size(sinogram);
  
 
    %% %%%%%%%%%%  initialize GPU %%%%%%%%%%%%%%%
    gpu  = gpuDevice();
    if ~isempty(par.GPU_list) && gpu.Index ~= par.GPU_list(1)
        % switch and !! reset !! GPU 
        gpu  = gpuDevice(par.GPU_list(1));
    end

     % ASTRA needs the reconstruction to be dividable by 32 othewise there
    % will be artefacts in left corner 
    Npix = ceil(Npix/r.binning);
    if isscalar(Npix)
        Npix = [Npix, Npix, Nlayers];
    elseif length(Npix) == 2
        Npix = [Npix, Nlayers];
    end
    
    if isempty(r.showed_layer_id)
        r.showed_layer_id = ceil(Npix(3)/2);
    end
    
    % !! important for binning => account for additional shift of the center
    % of rotation after binning, for binning == 1 the correction is zero
    rotation_center = [Nlayers, width_sinogram]/2; 

%     rotation_center(2) = rotation_center(2) + 0.5*(1-1/r.binning) ;
  
    if ~isempty(r.CoR_offset)
        rotation_center(2) = rotation_center(2) + r.CoR_offset/r.binning; 
    end
    
    if ~isempty(r.CoR_offset_v)
        rotation_center(1) = rotation_center(1) + r.CoR_offset_v/r.binning; 
    end
    
    % !! important for binning => account for additional shift of the center
    % of rotation after binning, for binning == 1 the correction is zero
    if par.is_laminography
        padding = 'symmetric'; 
    else  % Im really not sure why it differs from normal tomo, but I have it empirically tested 
        padding = 0; 
    end


    CoR_offsets_x = 0 ; 
    CoR_offsets_y = 0 ; 
    
    lamino_angles_offsets = 0 ; 
    tilt_angle_offsets = 0 ; 
    rot_angle_offsets = 0; 
    shear_angle_offsets = 0; 
    
    search_grid = linspace(r.search_range(1),r.search_range(2),r.num_grid_points); 
    switch lower(r.search_parameter)
        case 'center_of_rotation'
            CoR_offsets_x = search_grid; 
        case 'center_of_rotation_y'
            CoR_offsets_y = search_grid; 
        case 'lamino_angle'
            lamino_angles_offsets = search_grid; 
        case 'tilt_angle'
            tilt_angle_offsets = search_grid; 
        case 'rot_angle'
            rot_angle_offsets = search_grid; 
       case 'shear_angle'
            shear_angle_offsets = search_grid; 
        otherwise
            error('Missing option, choose from: center_of_rotation, lamino_angle, tilt_angle, rot_angle')
    end
    
    if par.usecircle && Npix(1) == Npix(2)
        radial_smooth_apodize= 10; 
        apodize = 20; 
        [~,circulo] = apply_3D_apodization(ones(Npix(1:2)), apodize, 0, radial_smooth_apodize); 
    end

    
    % generate dummy config 
     [cfg, vectors] = ...
            astra.ASTRA_initialize(Npix, [Nlayers, width_sinogram],angles ); 
    % use FBP function to provide already filtered sinogram 
    utils.verbose(0,'Filtering sinogram')
    [~,sinogram_filtered]  = FBP(sinogram, cfg, vectors, 1,...
        'GPU', par.GPU_list, 'verbose', 0, 'keep_on_GPU', true, ...
        'filter', par.filter_type, 'filter_value', par.freq_scale, ...
         'padding', padding, 'only_filter_sinogram', true);
    clear sinogram
    
%    plotting.smart_figure(1)
    clf     
    utils.verbose(0,'Parameter scan ... ')

    for ii = 1:length(search_grid)
        
        [cfg, vectors] = ...
            astra.ASTRA_initialize(Npix, [Nlayers, width_sinogram],...
                  angles        + r.rotation_angle_offset+rot_angle_offsets(min(ii,end)), ...
                r.lamino_angle_offset + par.lamino_angle + lamino_angles_offsets(min(ii,end)),...
                r.tilt_angle_offset + par.tilt_angle + tilt_angle_offsets(min(ii,end)), 1, ...
                 rotation_center + [(CoR_offsets_y(min(ii,end)))/r.binning,(CoR_offsets_x(min(ii,end)))/r.binning], ...
                r.shear_angle_offset + par.skewness_angle + shear_angle_offsets(min(ii,end)) ); 
                                                
        % find optimal split of the dataset for given GPU 
        split = astra.ASTRA_find_optimal_split(cfg, length(par.GPU_list),1,'back');

        %% backproject the already filtered sinogram method
        verbose(2,'FBP')
        rec  = tomo.Atx_sup_partial(sinogram_filtered, cfg, vectors, [1,1,length(par.GPU_list)],...
            'GPU', par.GPU_list, 'verbose', 0, 'split_sub', split);
        
        if par.usecircle && Npix(1) == Npix(2)
            rec = rec  .* circulo; 
        end    

        
        plotting.smart_figure(144)
        plotting.imagesc3D(rec, 'init_frame', r.showed_layer_id)
        axis image off
        colormap bone 
        title(sprintf('Lamino global param search: step id %i/%i', ii, length(search_grid)))
        drawnow 
        
        rec_preview_all(:,:,ii) = rec(:,:,max(1, min(end, r.showed_layer_id))); 
        
        
        [dX, dY] = math.get_img_grad(rec); 
        % estimate total variation 
        TV(ii) = gather(mean(mean2(abs(dX) + abs(dY)))); 
        STD(ii) = gather(std(rec(:))); 
        SP(ii) = gather(sparseness(abs(dX) + abs(dY)));
        utils.progressbar(ii, r.num_grid_points)
                
    end
    
    Nfine = 1e3; 
    fine_offsets = linspace(r.search_range(1),r.search_range(2),Nfine);
    
    TV = (TV - mean(TV)) / std(TV); 
    STD = (STD - mean(STD)) / std(STD); 
    SP = (SP - mean(SP)) / std(SP); 

    spline_TV = interp1(search_grid, TV, fine_offsets, 'spline'); 
    spline_STD = interp1(search_grid, STD, fine_offsets, 'spline'); 
    spline_SP = interp1(search_grid, SP, fine_offsets, 'spline'); 

    
    figure()
    hold all 
    plot(fine_offsets, spline_TV, '-r')
    plot(fine_offsets, spline_STD, '-b')
    plot(fine_offsets, spline_SP, '-G')
    plot(search_grid, TV, 'or')
    plot(search_grid, STD, 'ob')
    plot(search_grid, SP, 'oG')
    

    hold off 
    xlabel(['Required additional correction of ', r.search_parameter], 'Interpreter', 'none')
    ylabel('Value')
    legend({'Total variation', 'Standard deviation', 'Sparsity'})
    grid on 
    title(sprintf('Final score for global parameter search: %s', r.search_parameter), 'interpreter', 'none')
    drawnow 
    
    figure  
    plotting.imagesc3D(rec_preview_all)
    axis image off
    colormap bone 
    title(sprintf('Preview of reconstruction for all param steps: step id %i/%i', ii, length(search_grid)))
    drawnow 
        
  
end
% 
% function sinogram = unwrap_data(sinogram, method, boundary)
%     switch lower(method)
%         case 'fft_1d'
%             % unwrap the data by fft along slices 
%             sinogram = -math.unwrap2D_fft(sinogram, 2, boundary);
%         % case 'fft_2d'
%         %   % unwrap the data by 2D fft along slices 
%         %     sinogram = -math.unwrap2D_fft_split(sinogram, boundary);
%         case {'none', 'diff'}
% 
%         otherwise
%             error('Missing method')
%     end
% end

function spars = sparseness(x)
    %Hoyer's measure of sparsity for a vector
    % from  scipy.linalg import norm

    order_1 = 1;
    order_2 = 2;
    x = x(:);
    sqrt_n = sqrt(length(x));
    spars = (sqrt_n - norm(x, order_1) / norm(x, order_2)) / (sqrt_n - order_1);
end
    
function img = imreduce(img, ROI, binning)
    import math.*
    import utils.*
    
    isReal = isreal(img); 
    
    % crop the FOV after shift and before "binning"
    if ~isempty(ROI)
       
       img = img(ROI{:},:);   % crop to smaller ROI if provided  
                              % apply crop after imshift_fft
    end
    
    Np = size(img); 
    % perform FT interpolation instead of binning 
    img = interpolateFT_centered(img, ceil(Np(1:2)/binning/2)*2, -1);
    if isReal; img = real(img); end
end

