% weight_sino = estimate_reliability_region_grad(complex_projection, probe_size, subsample)
% Estimates the good region in the (real) sinogram for alignment and reconstruction.
% The method is base on gradient magnitude of the image 
% Written by YJ

function weight_sino = estimate_reliability_region_grad(sinogram, fill, erode_mat)
    %{
    disp(size(sinogram))
    [Gmag,~] = imgradient(gather(sinogram));
    sinogram_temp = imfill(Gmag,fill);
    level = graythresh(sinogram_temp);
    weight_sino_temp = imbinarize(sinogram_temp,level);
    weight_sino = imerode(weight_sino_temp,erode_mat);
    %}
    
    weight_sino = single(zeros(size(sinogram)));
    for i=1:size(sinogram,3)
        if ~isreal(sinogram)
            [Gmag,~] = imgradient(angle(sinogram(:,:,i)));
        else
            [Gmag,~] = imgradient(sinogram(:,:,i));
        end
        sinogram_temp = imfill(Gmag,fill);
        level = graythresh(gather(sinogram_temp));
        weight_sino_temp = imbinarize(gather(sinogram_temp),level);
        weight_sino(:,:,i) = imerode(weight_sino_temp,erode_mat);
    end
    
end
