% weights = estimate_reliability_region_grad(input, r_close, r_erode, par)
% Estimates good regions in real/complex projection stack base on gradient
% magnitude, with tricks such as morphological transformations.
% Inputs: 
%   **input       - real or complex image stack 
%   **r_close     - disk radius of structuring element for close transform
%   **r_erode     - disk radius of structuring element for erode transform
%   **par         - parameters to control other image processing steps
% *returns* 
%   ++output      - weights for each projection
% Written by YJ

function weights = estimate_reliability_region_grad(input, r_close, r_erode, par)

    weights = single(zeros(size(input)));
    
    r_close = double(gather(r_close));
    r_erode = double(gather(r_erode));
    
    for i=1:size(input,3)
        if ~isreal(input)
            sinogram_temp = angle(input(:,:,i));
        else
            sinogram_temp = input(:,:,i);
        end
        
        if par.unsharp
            sinogram_temp = imfilter(sinogram_temp, fspecial('unsharp'));
        end
        [sinogram_temp,~] = imgradient(sinogram_temp);
        if par.fill_connectivity>0
            sinogram_temp = imfill(sinogram_temp, par.fill_connectivity);
        end
        
        level = graythresh(gather(sinogram_temp));
        weight_temp = imbinarize(gather(sinogram_temp), level);
        if r_close(i)>0
            weight_temp = imclose(weight_temp,  strel('disk', r_close(i)));
        end
        if r_erode(i)>0
            weight_temp = imerode(weight_temp, strel('disk', r_erode(i)));
        end
        weights(:,:,i) = weight_temp;
    end
end
