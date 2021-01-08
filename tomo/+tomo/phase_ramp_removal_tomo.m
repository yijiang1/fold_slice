%  PHASE_RAMP_REMOVAL_TOMO Use tomography consistency between measured and 
%  reconstructed phase to remove phase ramp from data 
%  This function uses mask in volume space to accuratelly find regions of
%  air in the projection space. These regions are iterativelly forced
%  towards zero 
%
%  Several iterations are performed to further improve precision 
%
% object_full = phase_ramp_removal_tomo(object_full,object_ROI, theta, Npix, total_shift, par, varargin)
%
% Inputs:
%     **object_full         - complex-valued projections 
%     **object_ROI          - reliable region used for reconstruction 
%     **theta               - tomography angles 
%     **Npix                - size of reconstruction 
%     **par                 - tomography parameter structure 
%   *optional*  (or use values from par structure as default if provided)
%     **binning = 4                         - bin data to make reconstruction faster & more robust 
%     **positivity = true                   - apply positivity constaint
%     **auto_weighting = true               - give less weight to thic regions of the sample 
%     **fourier_guess = true                - calculate FFT to find phase ramp, important if the phase ramp is more than 2pi per frame 
%     **Niter = 3                           - number of iterations for phase removal
%     **unwrap_data_method = 'fft_2d'       - fft_2d , fft_1d
%     **sino_weights = 1                    - importance weights 
%     **CoR_offset = []                     - offset of the center of rotation, default is center of projection 
%     **inplace_processing = false          - process data inplace to save memory
%
% *returns* 
%    ++object_full          - complex-valued projections after phase ramp removal
    
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

function [object_full, W] = phase_ramp_removal_tomo(object_full,object_ROI, theta, Npix,total_shift, par, varargin)
    
    import utils.*
    verbose(struct('prefix', 'phase_ramp_remove'))

    parser = inputParser;
    parser.addParameter('binning', 4 , @isnumeric )     % bin data to make reconstruction faster & more robust 
    parser.addParameter('positivity', true , @islogical )  % apply positivity constaint
    parser.addParameter('auto_weighting', true , @islogical )  % give less weight to thic regions of the sample 
    parser.addParameter('fourier_guess', true , @islogical )  % calculate FFT to find phase ramp, important if the phase ramp is more than 2pi per frame 
    parser.addParameter('Niter', 3 , @isnumeric )  % number of iterations for phase removal
    parser.addParameter('unwrap_data_method', 'fft_2d' , @isstr )  % fft_2d , fft_1d
    parser.addParameter('sino_weights', 1, @isnumeric )  % importance weights 
    parser.addParameter('CoR_offset', [] , @isnumeric ) % offset of the center of rotation, default is center of projection 
    parser.addParameter('CoR_offset_v', [] , @isnumeric ) % added by YJ. vertical offset of the center of rotation, default is center of projection 
    parser.addParameter('inplace_processing', false, @islogical )  % process data inplace to save memory

    parser.parse(varargin{:})
    r = parser.Results;

    % load all varargins to the param structure 
    for name = fieldnames(r)'
        if ~isfield(par, name{1}) || ~ismember(name, parser.UsingDefaults) % prefer values in param structure if parsers returns default value 
            par.(name{1}) = r.(name{1});
        end
    end
    
    verbose(0,'Calculating phase ramp + amplitude correction')

   
     
    Np_full = size(object_full); 
 
    
    verbose(1,'Binning: %i', par.binning)
    
    if ismember(lower(par.unwrap_data_method), {'none', 'fft_2d'})
        interp_sign = -1 ;
    else
        interp_sign = 1 ;
    end
    
    % use symmetrically expanded ROI, get more region around sample 
    for ii = 1:2
        object_ROI{ii} = max(1, object_ROI{ii}(1)-ceil(par.asize(ii)/4)):min(Np_full(ii), ceil(object_ROI{ii}(end)+par.asize(ii)/4));
    end
    
    % shift the projections back to the "after loading" positions -> avoid boundary problems when
    % the phase ramp removal is applied 
    % !! high accuracy downsampling and shift is not needed in this function !!
    object = tomo.block_fun(@imshift_generic,object_full, -total_shift, [], [], 1, object_ROI, par.binning, 'fft', interp_sign,struct('use_fp16', false));
    
    Npix = ceil(Npix / par.binning);

    if isscalar(Npix) 
       Nlayers = size(object,1);
       Npix = [Npix,Npix,Nlayers]; 
    end
    
    Ngpu = max(1,length(par.GPU_list));
 
    if ~isscalar(par.sino_weights) && ~isempty(par.sino_weights)
        sino_weights = tomo.block_fun(@imshift_generic,par.sino_weights, -total_shift, Np_full(1:2), [], 1, object_ROI, par.binning, 'linear', ...
        struct('use_GPU', true, 'full_block_size', Np_full));
    else
        sino_weights = 1;
    end
    if all(mean(mean(abs(sino_weights-mean(mean(sino_weights))))) < 1e-2)
        sino_weights = 1;
    else
        sino_weights = real(sino_weights ./ max(max(sino_weights)));
    end
%     if ismatrix(par.illum_sum)
%         sino_weights = sino_weights .* imshift_generic(par.illum_sum,[0,0],Np_full(1:2), [], 1, object_ROI, par.binning, 'linear');
%     end
    
    gamma_tot = 1; 
    gamma_x_tot = 0; 
    gamma_y_tot = 0; 
 
    [~,circulo] = apply_3D_apodization(ones(Npix), 0, 0, 10); 


for ii = 1:par.Niter 
    progressbar(ii, par.Niter)


    phase = tomo.block_fun(@unwrap_object,object,sino_weights, par, struct('use_fp16', false, 'verbose_level', 0)); 
    if par.positivity
        % "positivity" constraint, useful for normal tomo but it has to be false for laminography
        phase = min(0, phase);
    end
    
    [Nlayers,width_sinogram,~]=size(phase);
    
    % find rotation center so that it stays consistent after binning
    
    par.rotation_center = [Nlayers, width_sinogram]/2; 

    if ~isempty(par.CoR_offset)  % important for laminography 
        par.rotation_center(2) = par.rotation_center(2) + par.CoR_offset/par.binning; 
    end
    %added by YJ
    if ~isempty(par.CoR_offset_v)  % important for laminography 
        par.rotation_center(1) = par.rotation_center(1) + par.CoR_offset_v/par.binning; 
    end
    
    par.rotation_center = par.rotation_center - total_shift(:,[2,1])/par.binning; 
     
    [cfg, vectors] = astra.ASTRA_initialize(Npix,[Nlayers,width_sinogram],theta,par.lamino_angle,par.tilt_angle, [par.horizontal_scale ; par.vertical_scale]', par.rotation_center); 
    split = astra.ASTRA_find_optimal_split(cfg,Ngpu,1,'back');

    % get FBP reconstruction from the initial guess 

    rec  = -tomo.FBP(phase, cfg, vectors, [1,1,Ngpu], 'GPU', par.GPU_list, 'split_sub', split, 'verbose',0);
    clear phase
    
    rec = rec .* circulo; % remove effect of unmeasured regions around sample 
    
    if par.positivity
        % positivity constraint 
        rec = max(0, rec);
    end
    % find model projections for given reconstruction 
    split = astra.ASTRA_find_optimal_split(cfg,Ngpu,1,'fwd');

    proj  = tomo.Ax_sup_partial(rec, cfg, vectors, [1,1,Ngpu], 'GPU', par.GPU_list, 'split_sub', split ,'verbose',0);
        
    if par.auto_weighting
        %% zero weights to regions with sample compared to air regions 
        Thresh = graythresh(rec(:)); 
        % find roughly region where is only air 
        mask = single(rec < Thresh); 
        % find the corresponding region in the projection space 
        proj_mask  = tomo.Ax_sup_partial(mask, cfg, vectors, [1,1,Ngpu], 'GPU', par.GPU_list, 'split_sub', split ,'verbose',0);
                  
        proj_blank  = astra.Ax_partial(ones(Npix,'single'), cfg, vectors, [1,1,Ngpu], 'GPU', par.GPU_list, 'split_sub', split ,'verbose',0);
        
        % define corresponding mask 
        W =  ((abs(proj_mask - proj_blank) ./ proj_blank) < 1e-2) .* sino_weights; 
        %size(proj_blank)
        % try to estimate weights direclty from the projections -> just to account for
        %  case when mask == 0 everywhere
        W = W + 1e-1*exp(-abs(proj).^2 / mean(abs(proj(:))).^2 );
        
    else
        W = sino_weights; 
    end
    
    W([1,end],:,:) = 0; % avoid boundary effects 

    % find phase ramp so that the masked regions are zero, if not possible, just enforce
    % consistency between the object and projection 

    [object, gamma, gamma_x, gamma_y] = stabilize_phase(object, exp(-1i*proj.* (1-W)), W, 'fourier_guess', false);

    
    gamma_tot = gamma_tot .* gamma; 
    gamma_x_tot = gamma_x_tot + gamma_x; 
    gamma_y_tot = gamma_y_tot + gamma_y;

end

if any(isnan(gamma_tot)) || any(isnan(gamma_x_tot)) || any(isnan(gamma_y_tot))
    error('Phase removal would result in NaNs')
end
if par.auto_weighting
    % store the produced mask -> false for regions of air 
    projection_mask =  ((abs(proj_mask - proj_blank) ./ proj_blank) < 1e-2) & (sino_weights > 0); 
end

%% calculate amplitude correction 
%use median of the masked regions to estimate amplitude correction factor
% use of median means the mask needs to be correct only in > 50% of the area
aobject = abs(object); 
if par.auto_weighting
    aobject(~projection_mask) = nan; 
end
amp_correction = reshape(nanmedian(reshape(aobject,[],Np_full(3))),1,1,[]);
% just to be sure that there is some mask everywhere
amp_correction(isnan(amp_correction)) = mean(mean(abs(object(:,:,isnan(amp_correction)))));


%% apply the phase and amplitude correction to the original stack_object array 
verbose(0,'Applying phase ramp + amplitude correction')

% Run locally on CPU , too slow GPU upload  / download 
cfg = struct('verbose_level',1,'inplace', par.inplace_processing, 'use_GPU', true); 
object_full = tomo.block_fun(@apply_ramp_shifted,object_full, gather(gamma_tot), gather(gamma_x_tot)/par.binning, gather(gamma_y_tot)/par.binning,total_shift,amp_correction, cfg);
verbose(0,'Done')
verbose(struct('prefix', 'template'))

    
end


%%% AUXILIARY FUNCTION FOR FAST PROCESSING ON GPU

function phase = unwrap_object(object,sino_weights, par)
    % get initial guess 
    switch lower(par.unwrap_data_method)
        case 'none' %added by YJ
            phase = angle(object);
        case 'fft_1d'
            phase = math.unwrap2D_fft(object,2,par.air_gap/par.binning);
        case 'fft_2d'
            phase = math.unwrap2D_fft2(object,par.air_gap/par.binning,0,sino_weights,1);
        otherwise
            error('Undefined unwrapping method')
    end
end


%%%% AUXILIARY FUNCTION FOR PARALLEL GPU PROCESSING  
               
function object_full = apply_ramp_shifted(object_full,gamma, gamma_x, gamma_y, total_shift, amp_correction)
    %  shift the projection to the original (ie after loading) positions to around ramp artefacts
    %  around edges if the projection was shifted too much 
    

    % it needs 2D circular shift (is nearest neighbor interpolation), FFT is not needed 
    object_full = utils.imshift_linear(object_full, -total_shift(:,1),-total_shift(:,2), 'circ'); 
  
    
    [M,N,~] = size(object_full);
    xramp = pi*(linspace(-1,1,M))';
    yramp = pi*(linspace(-1,1,N));
    if ~isa(object_full, 'gpuArray')
         object_full = auxfun(object_full, gamma, gamma_x, gamma_y, xramp, yramp, amp_correction);
    else
        % use inplace GPU calculation
        object_full = arrayfun(@auxfun, object_full, gamma, M*gamma_x, N*gamma_y, xramp, yramp, amp_correction);
    end
    
    object_full = utils.imshift_linear(object_full, total_shift(:,1),total_shift(:,2), 'circ');

end

function  object = auxfun(object, gamma, gamma_x, gamma_y, xramp, yramp, amp_correction)

    object = object .* (gamma./ amp_correction);  % correct global phase and also amplitude
    object = object .* exp(1i*xramp.*gamma_x); 
    object = object .* exp(1i*yramp.*gamma_y); 

end


