% simple and fast phase ramp/offset removal from already unwrapped sinogram 
% Assume the sample to contain air on both horizontal sides in stripes: 
% 1:air_gap(1) and end-air_gap(2):end
%
% sinogram = remove_sinogram_ramp(sinogram,air_gap,polyfit_order )
% 
% using linear interpolation between air on both sides 
% Inputs:
%     **sinogram - unwrapped projections 
%     **air_gap = [2x1] vector with number of pixels on both sides where can
%               be assumed to be air 
%     **polyfit_order:-1 = dont assume anything about the removed phase,
%                       subtract linear offset from each horizontal line separately
%                      0 = assume that removed degree of freedom is only a constant offset
%                      1 = assume that removed degree of freedom is a 2D plane 
% Outputs: 
%     ++unwrapped phase after the ramp removal 


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


function sinogram = remove_sinogram_ramp(sinogram,air_gap, polyfit_order)
    if nargin < 3
        polyfit_order = -1; 
    end 
    
    air_gap  = ceil(air_gap);
    [Nlayers,width_sinogram,~] = size(sinogram);
    ax = 1:width_sinogram;
    mask{1} = ax <= air_gap(1);
    mask{2} = ax >= width_sinogram-air_gap(min(end,2));
    %%%%%%%% get average sinogram to be subtracted 
    for ii = 1:2
        % get average values in the air_gap
        air_values{ii}  = sum( bsxfun(@times,sinogram,mask{ii}),2)/sum(mask{ii});
        if polyfit_order == 0
            %% find phase offset value 
            air_values{ii} = mean(air_values{ii}); 
        elseif polyfit_order == 1
            %% find 2D plane that fits the air region
            ramp = linspace(-1,1,Nlayers)';
            weight = 1; 
            % iterativelly refine the ideal plane estiamte to ignore
            % imperfect estimations of the air_gap paramter 
            for jj = 1:10
                % fit it with a linear plane, use for 2d unwrapping
                plane_fit = mean(weight .*air_values{ii}) ./ mean(weight) +mean(weight.*air_values{ii}  .* ramp)./ mean(weight.*ramp.^2)  .*ramp; 
                deviation = 5*mad(air_values{ii}(:) - plane_fit(:),0); 
                % avoid outliers by adaptive weighting
                weight = 1./(1+(abs(air_values{ii}-plane_fit)./deviation).^2  ); 
            end
            air_values{ii}  = plane_fit; 
        end
    end
    % get a plane in between left and right side passing through
    % air_values{1} and air_values{2}
    ramp = permute(interp1([0,width_sinogram]', permute([air_values{:}],[2,1,3]), 1:width_sinogram), [2,1,3]);
    %%%%%%%% remove ramp & sinogram offset
    sinogram =  sinogram - ramp; 


end