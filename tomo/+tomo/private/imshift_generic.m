
% IMSHIFT_GENERIC auxiliar function to be performed bu block_fun on GPU 
% it applies imshift_fft on the provided image that was first upsampled
% to Npix (if needed) and cropped to region ROI 
% after shifting, the image is downsampled by the chosen interpolation
% method "intep_method" that is more accurate than simple binning 
% 
% img = imshift_generic(img, shift, Npix, affine_matrix, smooth, ROI, downsample, intep_method, interp_sign)
% 
%  Inputs: 
%   **img           2D stacked image 
%   **shift         Nx2 vector of shifts applied on the image 
%   **Npix          2x1 int, size of the img to be upsampled before shift, Npix = [] -> no upsampling 
%   **affine_matrix affine metrix ! not implemented yet!
%   **smooth        how many pixels around edges will be smoothed before shifting the array 
%   **ROI           cell array, used to crop the array to smaller size 
%   **downsample    downsample factor , 1 == no downsampling
%   **intep_method  interpolation method: linear, fft 
%   **interp_sign   sign used for subpixel shifts of the dataset, +1 for unwrapped phase, -1 for phase differene 
%  *returns* 
%   ++img           2D stacked image 

function img = imshift_generic(img, shift, Npix, affine_matrix, smooth, ROI, downsample, intep_method, interp_sign)

    if nargin < 9 
        interp_sign = 0; 
    end
    
    import math.*
    import utils.*
    
    if isa(img, 'uint8') || (isa(img, 'gpuArray') && strcmpi(classUnderlying(img),'uint8'))
        img = single(img) / 255; % assume that the provided image is only compressed into uint8
    end
    
    % if needed upsample to the size of the projection 
    if ~isempty(Npix) && any(Npix(1:2) ~= [size(img,1),size(img,2)])
        switch intep_method
            case 'linear', img = utils.interpolate_linear(img, Npix); 
            case 'fft', img = utils.interpolateFT(img, Npix); 
        end
    end
    
    isReal = isreal(img); 
   
    if any(shift(:) ~=0 )
        smooth_axis = 3-find(any(shift ~= 0)); 
        img = smooth_edges(img, smooth, smooth_axis);
        if ~ismatrix(img)
            switch intep_method
                case 'linear'
                    img = utils.imshift_linear(img,shift);  % interpolation of the weights does not need such precision 
                case 'fft'
                    %%% APPLY SHIFT USING FFT -> periodic boundary 
                    img = imshift_fft(img, shift);            
            end
        end
    end

    % crop the FOV after shift and before "downsample"
    if ~isempty(ROI)
       img = img(ROI{:},:);   % crop to smaller ROI if provided  
                              % apply crop after imshift_fft
    end
    
    Np = size(img); 
    
    % perform interpolation instead of downsample , it provides more accurate results 
    if downsample > 1
        img = utils.imgaussfilt3_conv(img,[downsample,downsample,0]); 
        % correct for boundary effects of the convolution based smoothing
        img = img ./ utils.imgaussfilt3_conv(ones(Np(1:2), 'like', img),[downsample,downsample,0]); 
        switch intep_method
            case 'linear', img = utils.interpolate_linear(img,ceil(Np(1:2)/downsample/2)*2);  % interpolation of the weights does not need such precision 
            case 'fft',    img = utils.interpolateFT_centered(utils.smooth_edges(img, 2*downsample),ceil(Np(1:2)/downsample/2)*2, interp_sign); % accurate interpolation using FFT
        end
    end
    if isReal; img = real(img); end
    
end
