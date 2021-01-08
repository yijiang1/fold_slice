% PROJECTION_PROPAGATION_OPTIMIZATION Estimate optimal propagation distance to minimize amplitude
%
%  optimum = projection_propagation_optimization( stack_object, angles, range, ROI, par)
%
% Inputs:
%   **stack_object  - complex projections 
%   **angles        - projection angles 
%   **range         - scanning  range 
%   **ROI           - region of interest, cell 
%   **par           - parameter structure 

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


function optimum = projection_propagation_optimization( stack_object, angles, range, ROI, par)

    disp('Estimation of optimal focus')

    propagation_score = tomo.block_fun(@scan_propagation, stack_object, range, par, struct('use_fp16', false,'use_GPU', true, 'ROI', {ROI}, 'GPU_list', par.GPU_list));

    propagation_score = propagation_score - mean(mean(propagation_score,1),3); 
    propagation_score = propagation_score ./ std(std(propagation_score,[],1),[],3); 
    
    score = squeeze(trimmean(propagation_score,10,'round',3)); 
    
    subplot(1,3,1)
    plot(range'*1e6,squeeze(propagation_score(:,1,:)) , '-')
    title('Variance amplitude')
    xlabel('Propagation distance [\mum]')
    ylabel('Normalized local variance')
    grid on 
    
    hold all 
    plotting.vline(1e6*range(math.argmin(score(:,1))))
    hold off 
    subplot(1,3,2)
    plot(range*1e6,squeeze(propagation_score(:,2,:)) , '-')
    title('Variance phase')
    xlabel('Propagation distance [\mum]')
    ylabel('Normalized local variance')
    grid on 
    hold all 
    plotting.vline(1e6*range(math.argmax(score(:,2))), 'r:', 'Optimal propagation')
    hold off 
    
    
    optimum = sort([math.argmax(score(:,2)),math.argmin(score(:,1))]);
    optimum = 1e6*range(optimum);
    %suptitle(sprintf('Optimal propagation %3.1f - %3.1f um',optimum ))
    sprintf('Optimal propagation %3.1f - %3.1f um',optimum )
    propagation_score(:,2,:) = -propagation_score(:,2,:); 
    
    propagation_score = propagation_score ./ min(propagation_score,[],1); 

    [optim_shift_amp,ind] = find(squeeze(propagation_score(:,1,:)) == 1); 
    [~,uind] = unique(ind); 
    optim_shift_amp = optim_shift_amp(uind); 

    [optim_shift_phase,ind] = find(squeeze(propagation_score(:,2,:)) == 1); 
    [~,uind] = unique(ind); 
    optim_shift_phase = optim_shift_phase(uind); 

    subplot(2,3,3)
    plot(1e6*range(optim_shift_amp)+randn(size(optim_shift_amp))'*0.01, 1e6*range(optim_shift_phase)+randn(size(optim_shift_amp))'*0.01, 'o');
    title(sprintf('Correlation between phase/amplitude %3.2f', corr(optim_shift_amp, optim_shift_phase)))
    axis equal square 
    grid on 
    xlabel('Optimal shift from amplitude')
    ylabel('Optimal shift from phase')
    
    
    subplot(2,3,6)
    hold all 
    plot(angles, 1e6*range(optim_shift_amp), '.')
    plot(angles, 1e6*range(optim_shift_phase), '.')
    hold off 
    xlabel('Angles [deg]')
    ylabel('Optimal offset [\mum]')
    legend({'Amplitude', 'Phase'})
    axis tight  
    grid on 
    
%     optimum = (range(optim_shift_amp) + range(optim_shift_phase))/2; 
    optimum = range(optim_shift_amp); 

end


function variance =  scan_propagation(stack_object, range, par)

    Nproj = size(stack_object, 3); 
    for kk = 1:length(range)
        shift = range(kk);
        stack_object_prop = utils.prop_free_nf(stack_object, par.lambda, shift, par.pixel_size); 
        stack_object_amp   = abs(stack_object_prop); 
        stack_object_phase = -math.unwrap2D_fft2(stack_object_prop,par.air_gap,0,1,0);
        clear stack_object_prop
        
        
        % estimate local variance for amplitude 
        stack_object_amp = stack_object_amp-utils.imgaussfilt2_fft(stack_object_amp,3); 
        stack_object_amp = reshape(stack_object_amp,[],Nproj); 
        variance(kk,1,:) = std(stack_object_amp);
        
        % estimate local variance for phase
        stack_object_phase = stack_object_phase-utils.imgaussfilt2_fft(stack_object_phase,3); 
        stack_object_phase = reshape(stack_object_phase,[],Nproj); 
        variance(kk,2,:) = std(stack_object_phase);
    end
    
end
