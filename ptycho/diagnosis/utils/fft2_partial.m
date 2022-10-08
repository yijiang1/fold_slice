% FFT2_PARTIAL apply fftn only on smaller blocks (important for GPU)
% 
% x = fft2_partial(x,split)
% 
% Inputs: 
%   **x     - input array 
%  *optional*
%   **split - number of blocks to split the array before FFT to save the memory 
% 
% *returns*
%  ++x       fft2 transformed array 

    
    
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



function x = fft2_partial(x,split, inverse)
    if nargin < 3
        inverse = false; 
    end
    if nargin < 2 || isempty(split)
        %% empirical condition assuming FFT involved, may be too pesimistic 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mem_req =  numel(x)*8 * log2(max(size(x,1),size(x,2))); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isa(x,'gpuArray')
            x = complex(x);  % move directly to complex to include the expected memore requirements 
            gpu = gpuDevice; 
            split = ceil(mem_req / gpu.AvailableMemory); 
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
    
    
     if any(split > 1)
         Ndims = ndims(x);
         for ax = 1:2
            x = math.fft_partial(x,ax,1+mod(ax+1,Ndims), split, inverse);
         end  
     else
        if ~inverse
            x = fft2(x);
        else
            x = ifft2(x);
        end
     end
end
