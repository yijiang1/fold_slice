% output = apply_reliability_region(input, weight)
% Apply weighting functions to a stack of projections

% Inputs: 
%   **input             - real or complex image stack 
%   **weight            - weights for each projection
% 
% *returns* 
%   ++output            - projections multiplied by their weights

% Written by YJ

function output = apply_reliability_region(input, weight)
    output = single(zeros(size(input)));
    for i=1:size(input,3)
        output(:,:,i) = gather(input(:,:,i)) .* gather(weight(:,:,i));
    end
end
