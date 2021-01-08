% APPLY_PROBE_CONTRAINTS apply support constrains on the probe in the real space, fourier space or any other
% plane if provided, it useds factor in mode.support_back_propagation_factor to perform ASM propagation 
%
% probe = apply_probe_contraints(probe, mode)
%
% ** probe     complex array with probe / probes 
% ** mode       structure containing parameters for selected probe mode 
%
% returns:
% ** probe     complex array with probe / probes 
%
%


function probe = apply_probe_contraints(probe, mode)
    import math.*
    import utils.*
    import engines.GPU.shared.*

    if  ~isempty(mode.probe_support)
        % apply support contraint in real space (ir nearfield propagated )
        if ~isempty(mode.support_fwd_propagation_factor)
            if isscalar(mode.support_fwd_propagation_factor) && isinf(mode.support_fwd_propagation_factor)
                probe =  fftshift_2D(fft2(fftshift_2D(probe)));  % propagate to infinity
            else
                probe =  ifft2(fft2(probe) .* mode.support_propagation_factor);     
            end
        end   

        %% apply real-space support 
        probe = probe .* mode.probe_support;

        if ~isempty(mode.support_back_propagation_factor)
            if isscalar(mode.support_back_propagation_factor) && isinf(mode.support_back_propagation_factor)
                probe =  fftshift_2D(ifft2(fftshift_2D(probe)));  % propagate to infinity
            else
                probe =  ifft2(fft2(probe).* mode.support_back_propagation_factor);  
             end
        end   
    
    end
    
    Np_p = size(probe); 

    if mode.probe_scale_upd(end) > 0  && ~isempty(mode.probe_scale_window)
        % apply windowing to avoid boundary issues when subpixel probe
        % rescaling is used 
        probe = probe .* mode.probe_scale_window ; 
    end
    if  ~isempty(mode.probe_support_fft) || mode.probe_scale_upd(end) ~= 0 
        %% apply contraint in the detector plane
        
        % propagate probe on the detector 
        probe = fwd_fourier_proj(probe, mode);  

        if ~isempty(mode.probe_support_fft)
            probe = probe .* mode.probe_support_fft;
        end
        if mode.probe_scale_upd(end) < 0 && ~isempty(mode.probe_scale_window)
            probe = probe .*  mode.probe_scale_window ; 
        end
        
        % propagate probe back to the sample plane 
        probe = back_fourier_proj(probe, mode);  
    end    
end

