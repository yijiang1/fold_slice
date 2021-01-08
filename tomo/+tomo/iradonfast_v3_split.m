% tomogram = iradonfast_v3_split(p, theta, varargin)
% FUNCTION  wrapper around iradonfast_v3 that automatically splits the data into
% smaller blocks to avoid too large memory allocation 
% Inputs / Outputs: same as for iradonfast_v3

    
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



function tomogram = iradonfast_v3_split(p, theta, varargin)

    
    [Nlayers,tomo_size,Nangles] = size(p) ; 
    block_size = ceil(2e9/(tomo_size*Nangles*4)); % split the task into ~2GB blocks 
    Nblocks = ceil(Nlayers / block_size ); 

    % preallocate array to store results 
    tomogram = zeros(tomo_size, tomo_size, Nlayers, 'single');
    
    for ii = 1:Nblocks
        ind = 1+(ii-1)*block_size:min(ii*block_size, Nlayers);
        tomogram(:,:,ind) =tomo.iradonfast_v3(permute(p(ind,:,:),[2,3,1]), theta, varargin{:});   % Calculate slice
    end
end
