% FFT_PARTIAL apply fft only on smaller blocks (important for GPU)
% 
%  x = fft_partial(x,fft_axis,split_axis, split, inverse = false)
% 
% Inputs: 
%   **x         - input array 
%   **fft_axis  - axis along which is performed FFT 
%   **split_axis - axis along which is the array split 
%
%  *optional*
%   **split - number of blocks to split the array before FFT to save the memory 
%   **inverse - if true, use ifft intead of fft 
% 
% *returns*
%   ++x       fft transformed array 

    
    
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


function x = fft_partial(x,fft_axis,split_axis, split, inverse)
    persistent current_gpu 
    
    import math.*
    
    if nargin < 5
        inverse = false; % do inverse fft 
    end
    
    if nargin < 4 || isempty(split)
        % auto estimation of the split 
        mem_req = numel(x)*8 * log2(size(x,fft_axis)); 
        if isa(x,'gpuArray')
            % move directly to complex to include the expected memore requirements 
            x = complex(x);
            if isempty(current_gpu) || isnan(current_gpu.AvailableMemory)
                current_gpu = gpuDevice; 
            end
            split = ceil(mem_req / current_gpu.AvailableMemory); 
        else
            if mem_req < 50e9
                % assume to have at least 50GB ram free ...
                split = 1; 
            else
                avail_mem = utils.check_available_memory * 1e6;
                split = ceil(mem_req /avail_mem); 
            end
        end
    end


    
    if all(split == 1)
        if inverse
           x = ifft(x,[],fft_axis);
        else
           x = fft(x,[],fft_axis);
        end
        return 
    end
    
    
    % check GPU memory 
    Np = size(x);
    Nps = Np;
    Nps(split_axis) = ceil(Nps(split_axis)  / split);

    ind = {':', ':',':'};
    for i = 1:split
        ind{split_axis} = (1+(i-1)*Nps(split_axis)): min(Np(split_axis), i*Nps(split_axis));
        x_tmp = x(ind{:});
        if ~inverse
            x_tmp = fft(x_tmp, [], fft_axis);  % avoid additonal memory assignment
        else
            x_tmp = ifft(x_tmp, [], fft_axis); 
        end
        x(ind{:}) = x_tmp;
    end
    
end
