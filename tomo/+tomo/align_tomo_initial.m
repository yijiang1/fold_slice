% ALIGN_TOMO_INITIAL  Get fast initial guess of the vertical and horizontal alignment 
%
% [optimal_shift] = align_tomo_initial(stack_object, shift_init, angles, param, varargin)
%
% Inputs:
%     **stack_object        - complex-valued input array that will be unwrapped and used for alignment 
%     **shift_init          - initial guess of the shifts 
%     **angles              - corresponding angles (used only for sorting the projections)
% *optional*:   % if not provided, value from param is used as default 
%     **air_gap - empty region around sample where phase = 0 is assumed 
%     **vert_range - vertical range used for alignment , try to avoid highly
%                   scattering / residual features 
%     **phase_jumps_threshold - threshold above which the phase difference
%                            is assumed to be wrong and masked out 
%     **alignment_invariant - choose: phase_2D, phase_1D, phase_derivative, goldstein
%     **use_vertical_xcorr_guess - if true, use crosscorrelation for initial guess 
%     **data_filter - high pass filter constant 0=none, 0.005-0.02 seems to be optimal  
%     OTHER INPUTS DESCRIBED IN CODE 
%
% *returns* 
%     ++optimal_shift       - (Nangles x 1 array) = vertical shift to be applied on the stack_object in order to minimize the vertical mass fluctuation 


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


function [optimal_shift] = align_tomo_initial(stack_object, shift_init, angles, ROI, param, varargin)

    if nargin < 3
        param = struct();
    end

    parser = inputParser;
    parser.addParameter('vert_range', [] , @isnumeric )     % vertical range used for alignment 
    parser.addParameter('air_gap',  [50, 50] , @isnumeric ) % rough estimate of air region 
    parser.addParameter('phase_jumps_threshold',  1 , @isnumeric )    % threshold above which the phase difference is assume to be wrong 
    parser.addParameter('alignment_invariant',  'phase_2D' , @isstr )  % name of the invariant used for alignment 
    parser.addParameter('showsorted', true , @islogical )  % if the projections should be plotted sorted by angle 
    parser.addParameter('use_vertical_xcorr_guess', true , @islogical )  % get an initial guess by Xcorr, -> avoid trapping in local minima 
    parser.addParameter('data_filter', 0.02 , @isnumeric )  % high pass filtering to remove low spatial freq. errors 

    %% internal variables, usually no need to change 
    parser.addParameter('outer_loops_refinement',  3 ,  @isnumeric )  % number of outer loops for linear refinement step 
    parser.addParameter('N_SVD_modes',  10 , @isnumeric )   % number of SVD modes used to fill empty gaps in phase invariant 
    parser.addParameter('weights',  [] , @isnumeric )   % numeric of logical array contaning weights for each projection and each pixels for 2D unwrapping 
    parser.addParameter('windowautopos',  true , @islogical )   % distribute the plots over screen autimatically 

    parser.parse(varargin{:})
    r = parser.Results;

    % load all to the param structure 
    for name = fieldnames(r)'
        if ~isfield(param, name{1})  % prefer values in param structure 
            param.(name{1}) = r.(name{1});
        end
    end

    import tomo.*
    import utils.*
    import math.*
    utils.verbose(struct('prefix', 'align'))

    if ~isempty(param.vert_range)
        ROI{1} = ROI{1}(max(1,param.vert_range(1)):min(end,param.vert_range(end)));
    end
    Nlayers = length(ROI{1}); 
    Nw = length(ROI{2}); 
    Nangles = length(angles);    

    switch param.alignment_invariant
        case 'phase_1D'
            %% standard vertical mass fluctuation
            utils.verbose(0,'Fast 1D FFT unwrapping')
            [phase, phase_diff, residues] = tomo.block_fun(@unwrap2D_fft, stack_object, 2, param.air_gap, struct('ROI', {ROI}, 'use_fp16', false));
            invar0 = max(0,squeeze(sum(-phase,2)));
            jumps = abs(phase_diff) > param.phase_jumps_threshold;
            mask_invar = squeeze(any(jumps,2));  %% relevance weights for each line 
            mask_residues = squeeze(sum(residues,2))>0; 
            mask_residues = conv2(mask_residues,ones(3,1),'same')>0;
            mask_invar(2:end,:) = mask_invar(2:end,:) | mask_residues; 
        case 'phase_2D'
            utils.verbose(0,'Fast 2D FFT unwrapping')
            % standard vertical mass fluctuation
            % fft-based unwrapping
            
            [phase, residues] = unwrap2D_fft2_split(stack_object, param.air_gap,1,param.weights,param.GPU_list,ROI);
            
            invar0 = squeeze(sum(phase,2));
            mask_invar = squeeze(sum(residues,2))>0; 
%             mask_invar = conv2(mask_invar,ones(3,1),'same')>0;
                  
        case 'phase_derivative'
            % vertical derivative fluctuation
            utils.verbose(0,'Get phase gradient')
            phase_diff = tomo.block_fun(@get_phase_gradient_1D,stack_object, 1,1, struct('ROI', {ROI}, 'use_fp16', false));
            invar0 = squeeze(sum(phase_diff,2));
            jumps = abs(phase_diff) > param.phase_jumps_threshold;
            jumps(:,[1:2,end-1:end],:) = 0; % avoid jumps caused by phase ramp 
            mask_invar = squeeze(sum(jumps,2) > 1);  %% relevance weights for each line 
        case 'phase_goldstein'
            utils.verbose(0,'Estimating residua')
            residues = abs(findresidues(stack_object)) > 0.1; 
            mask_invar = squeeze(sum(residues,2))>0; 
            if sum2(mask_invar) > 1
                warning('Selected range contains %i residua', sum2(mask_invar))
            end
            phase = zeros(Nlayers, Nw, Nangles, 'single'); 
            parfor ii = 1:Nangles
                utils.progressbar(ii, Nangles)
                o = stack_object(:,:,ii)
                phase(:,:,ii) = utils.goldsteinunwrap2(angle(o(ROI{:})));
            end
            phase = utils.remove_sinogram_ramp(phase,param.air_gap, true); 
            invar0 = squeeze(sum(phase,2));
        otherwise
            error('Missing option %s', par.alignment_invariant)
    end
    clear jumps 
    % move to GPU, always assume that GPU is availible 
    invar0 = Garray(invar0);
    
    if any(sum(invar0)==0)
       error('Some projections are empty') 
    end
    if ~exist('residues', 'var') && ~strcmpi(param.alignment_invariant, 'phase_derivative')
        utils.verbose(0,'Estimating residua')
        residues = tomo.block_fun(@(x)(abs(findresidues(x)) > 0.1), stack_object); 
        mask_invar = mask_invar | squeeze(sum(residues,2))>0; 
    end
    if any(mean(mask_invar) > 0.9)
       wrong = param.scanstomo(mean(mask_invar) > 0.9);
       error(sprintf(['Too many phase jumps in %i angles, alignment will fail \n try to increase par.phase_jumps_threshold or change par.alignment_invariant\n Wrong scans: ', repmat('%i ',1,length(wrong)), ' \n quitting'], length(wrong), wrong ))
    elseif any(mean(mask_invar) > 0.7)
       wrong = param.scanstomo(mean(mask_invar) > 0.7);
       warning(sprintf(['Too many phase jumps in %i angles, alignment will most likely fail \n try to increase par.phase_jumps_threshold or change par.alignment_invariant\n Wrong scans: ', repmat('%i ',1,length(wrong)), ], length(wrong), wrong ))
    end
    
    if isempty(shift_init)
       shift_init = zeros(Nangles, 1);
    end
    Nplots = 2+param.use_vertical_xcorr_guess;
    
    if param.showsorted
        [sangles,plot_sort] = sort(angles);
        x_axis = sangles;
        x_label = 'Angled [deg]';
    else
        plot_sort = 1:Nangles;
        x_axis = param.scanstomo;
        x_label = 'Scan number';
    end
    weight = ~mask_invar; 
   

        % remove linear offset -> prevents boundary problems 
        invar = remove_linear_ramp(invar0);

        % apply only integer shift 
        shift_Y = shift_init;
        invar = imshift_fft_ax(invar, shift_Y, 1);
        weight = imshift_linear_ax(weight, shift_Y, 1, 'nearest', 0);
        % select range without boundary issues 
        offset= max(abs(shift_Y));
        range = round(2+offset : Nlayers - offset-1); 
        if length(range) < 20; error('Too small range for vertical alignment'); end
        invar = invar(range,:); 
        weight = weight(range,:);
        Nlayers = length(range);
        
        % remove linear offset -> prevents boundary problems 
        invar = remove_linear_ramp(invar);
        

        %% plot initial alignment 
        fig_id = 5667;
        if param.windowautopos && ~ishandle(fig_id)  % autopositioning only if the figure does not exists yet 
            plotting.smart_figure(fig_id)
            set(gcf,'units','normalized','outerposition',[0.2 0.2 0.8 0.8])
        else
            plotting.smart_figure(fig_id)
        end
    
        ax(1)=subplot(Nplots,3,1);
        invar_tmp = invar;
        invar_tmp = imfilter_high_pass_1d(invar_tmp, 1, param.data_filter, Nlayers/2);
        imagesc(x_axis, 1:Nlayers, invar_tmp(:,plot_sort), quantile(invar_tmp(~isnan(invar_tmp)), [1e-2,1-1e-2]))
        axis xy
        grid on 
        title('No alignment, linear ramp removed')
        ylabel('Vertical axis [pixels]')
        xlabel(x_label)
        subplot(Nplots,3,2)
        invar_tmp(~weight) = nan ;
        plot(invar_tmp)
        xlabel('Vertical pixels')
        axis tight 
        ax(4)=subplot(Nplots,3,3);
        imagesc(x_axis, 1:Nlayers, 1-weight(:,plot_sort))
        axis xy
        title('Phase jumps / Residues mask')
        utils.verbose(0,'Vertical alignment - initial guess')
        utils.verbose(0,'Damaged pixels: %3.2g%%', mean2(~weight)*100)
        %subtitle('Tomography invariant vertical alignment')
        ylabel('Vertical axis [pixels]')
        xlabel(x_label)
        
        if param.use_vertical_xcorr_guess
            %% use cross correlation as the first guess 

    
            [shift_Y,invar_filtered] = ...
                cross_correlation_estimation(invar, weight, angles, param.N_SVD_modes, param.data_filter);
                        
            % try to be smart and avoid drastic jumps 
%             shift_Y = max(shift_Y, quantile(shift_Y, 1e-2));
%             shift_Y = min(shift_Y, quantile(shift_Y, 1-1e-2));
            % minimize the shift offset 
            shift_Y = shift_Y - (max(shift_Y)+min(shift_Y))/2;
            %shift_Y = shift_Y - median(shift_Y);
                    
            % perform only nearest neighbor shift 
            invar_filtered = imshift_linear_ax(invar_filtered, shift_Y, 1, 'circ');
            weight_shifted = imshift_linear_ax(weight, shift_Y, 1, 'nearest',0);

            %% plot current estimation 
            ax(2)=subplot(Nplots,3,4);
            imagesc(x_axis, 1:Nlayers, invar_filtered(:,plot_sort), quantile(invar_filtered(:), [1e-2,1-1e-2]))
            title('X-corr based guess - highpass filtered')
            axis xy
            grid on 
            ylabel('Vertical axis [pixels]')
            xlabel(x_label)
            subplot(Nplots,3,5)
            invar_filtered(~weight_shifted) = nan;
            plot(invar_filtered)
            axis tight 
            xlabel('Vertical pixels')
            title('Line plot - highpass filtered')
            
            subplot(Nplots,3,6)
            plot(x_axis, shift_Y(plot_sort))
            axis tight 
            title('Applied shift')
            ylabel('Shift [pixels]')
            grid on 
            utils.verbose(0,'Vertical alignment - iterative refinement')
            xlabel(x_label)
        else
            shift_Y = zeros(Nangles,1);
        end
                
        %% iterative vertical position refinement 
        for ii = 1:param.outer_loops_refinement
            progressbar(ii, param.outer_loops_refinement)
            % shift sinograms 
            invar_shifted =  imshift_fft_ax(invar, squeeze(shift_Y),1);
            weights_shifted =  imshift_linear_ax(weight, squeeze(shift_Y),1,'nearest',0);

            % select range without boundary issues 
            offset= max(abs(shift_Y));
            range = round(1+offset : Nlayers - offset); 
            assert(length(range) > 30, 'Too small vertical range for alignment')
            % crop to the undamaged region by boundary issues 

            invar_shifted = invar_shifted(range,:); 
            weights_shifted = weights_shifted(range,:);
            % perform alignment 
            [shift_update, invar_shifted,weights_shifted] = linear_iterative_refinement(invar_shifted, weights_shifted, param.data_filter);
            shift_Y = shift_Y + shift_update;
        end
        
        if param.outer_loops_refinement > 0        
            %% plot results 
            ax(3)=subplot(Nplots,3,3*Nplots-2);
            imagesc(x_axis, range , invar_shifted(:,plot_sort), quantile(invar_shifted(:), [1e-2,1-1e-2]))
            axis xy
            grid on 
            title('Iterative refinement')
            xlabel(x_label)
            ylabel('Vertical axis [pixels]')

            %% plot results 
            subplot(Nplots,3,3*Nplots-1)
            invar_shifted_plot = invar_shifted;
            invar_shifted_plot(weights_shifted==0) = nan;
            plot(invar_shifted_plot)
            xlabel('Vertical pixels')
            axis tight 
            title('Line plot - highpass filtered')
            subplot(Nplots,3,3*Nplots)
            plot(x_axis, shift_Y(plot_sort))
            axis tight 
            title('Applied shift')
            ylabel('Shift [pixels]')
            grid on 
            xlabel(x_label)
        end
        
        try linkaxes(ax, 'xy'); end
        
        drawnow 
                
    if exist(param.output_folder, 'file') && ~debug()
        try
            paths{1} = [param.output_folder, '/vertical_alignment.png']; 
            if param.online_tomo
                paths{2} = [param.output_folder, '/vertical_alignment.png']; 
            end
            for ii = 1:length(ii)
                print(['-f', num2str(fig_id)],'-dpng', paths{ii} )
                utils.verbose(0,['Plot saved to:', paths{ii}])
                system(sprintf('convert -trim %s %s', paths{ii}, paths{ii}));
            end
        catch err
            warning('vertical_alignment.png saving failed with error:\n "%s"', err.message)
        end
    end
  
    
    optimal_shift = shift_Y + shift_init;
    
    optimal_shift = optimal_shift - median(optimal_shift);

end
function [shift_Y, invar, weight] = cross_correlation_estimation(invar0, weight, angles, N_SVD_modes, data_filter)
    %% cross-corelation based alignment guess 
    %% make a fast initial guess based on the tomography invariant and cross-correlation
    % take several angles at the beginning to get some initial guess of the
    % vertical fluctuation shape and use Xcorr to find the optimal shifts 
    
    import utils.Garray

    % remove linear offset 
    [Nlayers, Nangles]= size(invar0);
    [~,angle_sort] = sort(angles);
    
    % move on GPU 
    invar0 = Garray(invar0); 
    weight = Garray(weight); 
    
       
    % helps a lot in case of golden ratio datasets 
    invar0 = invar0(:,angle_sort);
    weight = weight(:,angle_sort);

    % apply high pass filter 
    invar = imfilter_high_pass_1d(invar0, 1, data_filter, Nlayers/2);
    % further suppress boundary effects 
    invar = invar.* tukeywin(Nlayers, 0.1);
        
    % update weight of inreliable pixels 
    range = quantile(invar(:), [0.01 , 0.99]); 
    weight_invar = invar > range(1) & invar < range(2) & weight;
   
    % crop to the limited range 
    invar = max(min(invar, range(2)), range(1));
    % fill missing data 
    invar = fill_gaps_1D(invar,~weight_invar, N_SVD_modes, 20);

    % find 5 of the most representative angles to be used as referene 
    [~, ~, ~, D] = kmeans(invar0',1);  % using invar before highpass filter to find optimal cluster center seems to work better 
    [~,ind] = sort(D);
    invar_reference = median(invar(:,ind(1:5)),2); 

    % find optimal shift using Xcorr method
    shift_Y = -utils.find_shift_fast_1D(invar,invar_reference,1,0)'; 

    
    % provide some extra robustness by using median filter -> assume that
    % neighboring projections are quite well aligned 
    shift_Y = gather(shift_Y);
    medfilt_win = 3; 
    mshift_Y = medfilt1(shift_Y,medfilt_win, 'truncate'); 
    medfilt_resid = shift_Y - mshift_Y; % residuum betw
    
    % avoid too large jumps with respect to rest of the shifts 
    range = 2*quantile(medfilt_resid, [0.001, 0.999]);
    medfilt_resid = max(min(medfilt_resid, range(2)), range(1));
    shift_Y = mshift_Y + medfilt_resid;
    
    weight = weight & weight_invar; 

    %% resort to original order 
    [~,scan_sort] = sort(angle_sort);
    shift_Y = shift_Y(scan_sort);
    weight = weight(:,scan_sort);
    invar = invar(:,scan_sort);
    
    % move from GPU 
    invar = gather(invar); 
    weight = gather(weight); 
    
    
end

function invar = fill_gaps_1D(invar0,mask_invar, N_SVD_modes, Niter)
     % try to repair failed values , iterativelly replace them using SVD
     % method by most propable value 
     
    import math.*
    if ~any(mask_invar(:))
        invar = invar0; 
        return; 
    end 
    invar = invar0; 
    range = quantile(invar(~mask_invar), [0.01, 0.99]); 

    for i = 1:Niter
        % slowly increase complexity 
        [U,S,V]=fsvd(invar,N_SVD_modes); 
        invar_filt = U * S*V';  % %get smooth estimate from SVD 
        invar = invar.*~mask_invar + invar_filt .* mask_invar ; %% replace missing by a smooth curve 
        % avoid outliers 
        invar = max(min(invar, range(2)), range(1)); 
    end 
    
end
function array = remove_linear_ramp(array)
    % auxiliary function to subtract linear ramp from sinogram 
    % it is important to avoid edge ringing and other artefacts when FFT
    % filtering is applied on the 2D array 
    
    [Nlayers]= size(array,1);
    Nedge = 5; % number of averaged  edge layers 
    top = mean(array(1:Nedge,:));
    bottom = mean(array(end-Nedge:end,:));
    ramp = interp1([0,Nlayers]',[top;bottom], 1:Nlayers);
    array =  array - ramp; 
end

function [total_shift_Y, invar,weights]  = linear_iterative_refinement(invar_0, weights_0, data_filter)
    %% ITERATIVE REFINEMENT OF VERTICAL ALIGNMENT METHOD 
    % method based on optical flow, it can deal better with the missing /
    % damaged data compared to the Xcorr based methods -> it us used for
    % refinement of the Xcorr guess 
    
    
    import utils.*
    import math.*

    [Nlayers,Nangles] = size(invar_0); 
        
    total_shift_Y = zeros(Nangles,1);

    invar_0 = Garray(invar_0); 
    weights_0 = Garray(single(weights_0)); 
    total_shift_Y = Garray(total_shift_Y);
        
    % apply high pass filter 
    invar_0 = remove_linear_ramp(invar_0);
    invar_0 = imfilter_high_pass_1d(invar_0, 1, data_filter, Nlayers/2);
    % further suppress boundary effects 
    invar_0 = invar_0 .* tukeywin(Nlayers, 0.1);
        
    % fill missing data 
    invar_0 = fill_gaps_1D(invar_0,~weights_0, 2, 20);


%      img = imshift_2D(ones(10), randn(10,1)*2)
     
     
    X = utils.Garray(1:Nlayers); 
    Y = utils.Garray(1:Nangles); 
    [X,Y] = meshgrid(X,Y);

     
    relax_step = 0.9; % avoid too large steps 
        
    for i = 1:1e3
        % run till convergence criterion is reached 
        
        invar =  imshift_fft_ax(invar_0, total_shift_Y, 1);  
        weights = interp2(weights_0, Y',X'+total_shift_Y', 'nearest', 0); 

        % apply high pass filter => get rid of phase artefacts 
        invar = imfilter_high_pass_1d(invar,1,data_filter, Nlayers/2); 
       
        % further suppress boundary effects 
        invar = invar .* tukeywin(Nlayers, 0.2);
        
        % take median over all the positions 
        m_invar = sum(invar .* weights,2) ./ (sum(weights,2)+1e-3); 
        % get gradient by convoolution to avoid edge issues when using fft
        md_invar = math.get_img_grad_conv(m_invar,2,1); 

        % in vertical direction use shift of the invariant => more robust and faster 
        DY = m_invar-invar; 
        
        shift_Y  = -squeeze(sum(weights .* (DY .* md_invar) ,1) ./...
                            sum(weights .* md_invar.^2,1));
        % avoid too large steps where linear approximation is not valid anymore
        shift_Y = relax_step * min(0.5,abs(shift_Y)) .* sign(shift_Y);
        total_shift_Y = total_shift_Y  + shift_Y'; 
                              
        err(i) = gather(mean2(weights .* DY.^2)); 
        if i > 2 && err(i-1) < err(i) || max(abs(shift_Y)) < 1e-2
            break   %  it will stop when alignment reaches the numerical precision 
        end
    end
    total_shift_Y = gather(total_shift_Y);
             
    invar = gather(invar);
    
            
end
