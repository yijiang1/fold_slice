% weight_sino = estimate_reliability_region_grad(complex_projection, probe_size, subsample)
% Estimates the good region in the (real) sinogram for alignment and reconstruction.
% The method is base on gradient magnitude of the image 
% Written by YJ

function weight_sino = estimate_reliability_region_grad(sinogram, fill, SE)

    weight_sino = single(zeros(size(sinogram)));
    
    for i=1:size(sinogram,3)
        if ~isreal(sinogram)
            sinogram_temp = angle(sinogram(:,:,i));
        else
            sinogram_temp = sinogram(:,:,i);
        end
        
        H = fspecial('unsharp');
        sinogram_temp = imfilter(sinogram_temp,H);
        [Gmag,~] = imgradient(sinogram_temp);
        sinogram_temp = imfill(Gmag,fill);
        level = graythresh(gather(sinogram_temp));
        
        weight_sino_temp = imbinarize(gather(sinogram_temp),level);
        weight_sino_temp = imclose(weight_sino_temp,SE);
        weight_sino_temp = imerode(weight_sino_temp,SE);
        weight_sino(:,:,i) = weight_sino_temp;
        
    end
end
