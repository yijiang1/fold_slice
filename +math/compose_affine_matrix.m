% COMPOSE_AFFINE_MATRIX calculate affine matrix when provided rotation, shear, asymmetry and scale 
% 
%  affine_mat  = compose_affine_matrix(scale, asymmetry, rotation, shear)
%  
% Inputs:
% **scale      A1 = [scale, 0; 0, scale]
% **asymmetry   A2 = [1+asymmetry/2,0; 0,1-asymmetry/2]
% **rotation    A3 = [cosd(rotation), sind(rotation); -sind(rotation), cosd(rotation)]
% **shear       A4 = [1,0;tand(shear),1];
% 
% returns: 
% ++ affine_mat  affine matrix = A1*A2*A3*A4

function affine_mat  = compose_affine_matrix(scale, asymmetry, rotation, shear)
    if isscalar(scale) && isscalar(asymmetry) && isscalar(rotation) && isscalar(shear)
        affine_mat = scale(1)*[1+asymmetry/2,0; 0,1-asymmetry/2]*[cosd(rotation), sind(rotation); -sind(rotation), cosd(rotation)] * [1,0;tand(shear),1];
    else
        for ii = 1:max([numel(scale), numel(asymmetry), numel(rotation), numel(shear)])
            affine_mat(:,:,ii) = scale(min(ii,end))*...
                                [1+asymmetry(min(ii,end))/2,0; 0,1-asymmetry(min(ii,end))/2]*...
                                [cosd(rotation(min(ii,end))), sind(rotation(min(ii,end))); -sind(rotation(min(ii,end))), cosd(rotation(min(ii,end)))] *...
                                [1,0;tand(shear(min(ii,end))),1];
        end
    end
end
