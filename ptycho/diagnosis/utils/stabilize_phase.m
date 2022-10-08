% [img_out, gamma,gamma_x, gamma_y, c_offset] = stabilize_phase(img, varargin)
% Description:  adjust phase of the object to be mostly around zero or close 
% to the provided reference image img_ref and remove linear ramp, phase
% offset and if normalize_amplitude==true normalize amplitude to be around 1
%
% Method: 
%  1) if fourier_guess == true and remove_ramp == true, get rough estiamte of the center of FFT and
%  use it to subtract phase ramp, IT CAN FAIL IF SAMPLE HAS STRONG
%  AMPLITUDE AND PHASE VARIATION, IT IS NOT ABLE TO USE THE WEIGHTS 
%  2) if remove_ramp == true, accuratelly refine the phase ramp by weighted
%  LSQ method, regions were W is small have small importance in the phase
%  ramp estiamtion. SECOND STEP ASSUMES THAT THE PHASE RAMP IN IMAGE IS
%  SMALLER THAN 2PI ACROSS THE IMAGE
%
% Inputs:
%     **img_orig        - complex image to stabilize
%     **img_ref         - complex images used as reference, if empty ones is used
% *optional*:
%     **weights         - (array) 0<x<1 numeric array denoting reliable unwrapping region 
%     **split           - (int) split FFTon smaller blocks on GPU 
%     **fourier_guess 	- (bool), use maximum in Fourier space to get initial guess of the phase ramp 
%     **remove_ramp     - (bool), if false, remove only phase offset (scalar)
%     **binning         - (int) bin the array to speed up phase ramp calculation and make it more robust
%     **normalize_amplitude - (bool) normalize amplitude to have average value around one
%     **split           - (scalar)  split volume, important for GPU when fft of large array is calculated 
%
% returns: 
%     ++img_out - phase ramp / offset subtracted complex-valued image 
%     ++gamma,gamma_x, gamma_y - (1,1,N arrays) offset, ramp horizontal / vertical
%     ++c_offset - 2D/3D array - either constant or phase ramp offset subtracted directly from the
%                               complex data as img .* exp(1i*c_offset) 
%
%    Example how to correct output image using gamma,gamma_x, gamma_y
%         [M0, N0] = size(img_orig);
%         c_offset = angle(gamma) + xramp.*gamma_x*M0 + yramp.*gamma_y*N0;
%         img_out = img_orig.*exp(1i*c_offset);

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



function [img_out, gamma,gamma_x, gamma_y, c_offset] = stabilize_phase(img_orig, varargin)

    %import math.*
    %import utils.*
    
    par = inputParser;
    par.addOptional('img_ref',  [] )
    par.addOptional('weights',  [] )
    par.addParameter('fourier_guess',  true , @islogical )   % use cross correlation as initial guess, avoids local minima but can be dangerous
    par.addParameter('remove_ramp',  true , @islogical )     % remove also ramp, not only offset 
    par.addParameter('binning', 0, @isint )                  % bin the array to speed up phase ramp calculation 
    par.addParameter('normalize_amplitude', false, @islogical )% normalize amplitude to have average value around one
    par.addParameter('split', 1, @isnumeric )                % split volume, important for GPU when fft of large array is calculated 

    par.parse(varargin{:})
    r = par.Results;
    img_ref = r.img_ref; 
    weights = r.weights; 
    binning = max(1, r.binning); 
    
    [M0,N0,~] = size(img_orig);
    if isreal(img_orig)
        return  % if the input is real, no phase shift is needed 
    end
    if isempty(img_ref)
        img_ref = 1;
    end
    if isempty(weights)
        weights = 1;
    end
    
    if isvector(img_orig)
        %% calculate the optimal phase shift 
        gamma = mean2(img_ref .* conj(img_orig));
        gamma = gamma./abs(gamma);
        if isnan(gamma); gamma = 1; end
        img_out = img_orig * gamma;
        return
    end
    
    img = img_orig; 
        
    if binning > 1
       % speed up calculation and make more robust by binning 
       img    = utils.binning_2D(img, binning);
       if ~isscalar(img_ref) || isempty(img_ref)
            img_ref= utils.binning_2D(img_ref, binning);
       end
       if ~isscalar(weights) || isempty(weights)
            weights= utils.binning_2D(weights, binning);
       end
    end
    
    %% calculate complex phase difference
    phase_diff = img_ref .* conj(img);

    [M,N,~] = size(img);
    xramp = pi*(linspace(-1,1,M))';
    yramp = pi*(linspace(-1,1,N));
    
    if r.fourier_guess && r.remove_ramp
        %% initial guess based on position of maximum in Fourier domain
        xcorrmat = ifftshift_2D(abs(ifft2_partial((phase_diff), r.split))).^2; 
        % center of mass seems to be more accurate if the phase FFT has bimodal distribution, now I
        % take center of mass of regions > 0.5 of xcorr maximum 
        [y,x] = center(max(0,xcorrmat-max(max(xcorrmat))/2), false);
        x=x-ceil(M/2)-1;
        y=y-ceil(N/2)-1;

        c_offset  =  xramp.*x + yramp.*y;
        phase_diff = phase_diff.*exp(1i*c_offset);
        
    end
        
    %% calculate the optimal phase shift 
    gamma = mean2(phase_diff .*weights) ./ mean2(weights);
    gamma = gamma./abs(gamma);
    if any(isnan(gamma)); gamma = 1; end

    if r.remove_ramp
        phase_diff = phase_diff .* conj(gamma);

        %% linear refinement 
        phase_diff= angle(phase_diff).*weights; % linearize the problem 
        gamma_x=mean2(phase_diff.*xramp) ./ mean2(weights.*abs(xramp).^2);
        gamma_y=mean2(phase_diff.*yramp) ./ mean2(weights.*abs(yramp).^2);
        if r.fourier_guess
            %% get total correction 
            gamma_x = gamma_x - x;
            gamma_y = gamma_y - y;
        end
        
        % export dimensionless 
        gamma_x = gamma_x / M; 
        gamma_y = gamma_y / N; 
        
        %% correct output  image 
        xramp = pi*(linspace(-1,1,M0))';
        yramp = pi*(linspace(-1,1,N0));

        if isa(img_orig, 'gpuArray')
            xramp = gpuArray(xramp); 
            yramp = gpuArray(yramp); 
            [img_out,c_offset] = arrayfun(@remove_phase,img_orig, xramp, yramp, gamma, gamma_x*M0/binning, gamma_y*N0/binning);
        else
            [img_out,c_offset] = remove_phase(img_orig, xramp, yramp, gamma, gamma_x*M0/binning, gamma_y*N0/binning);
        end

     else
         img_out = img_orig .* gamma;
         c_offset = angle(gamma); 
    end
     
    if r.normalize_amplitude
        mean_amplitude = mean2(weights .* img_out) ./ mean2(weights); 
        img_out = img_out ./ mean_amplitude;
    end

end

% auxiliar function for GPU processing (kernel merging)
function [img_out,c_offset] = remove_phase(img_orig, xramp, yramp, gamma, gamma_x, gamma_y)

    %% correct output  image 
    c_offset = angle(gamma) + xramp.*gamma_x + yramp.*gamma_y;
    img_out = img_orig.*exp(1i*c_offset);


end