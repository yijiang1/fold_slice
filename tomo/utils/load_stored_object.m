% LOAD_STORED_OBJECT load stored complex projections, ( + angles and par structure) saved by function utils.savefast_save 
%
%  [stack_object, theta, par] = load_stored_object(path)
%
% Inputs: 
%   **path        -path to the stored object  
% *returns*
%   ++stack_object measured projections 
%   ++theta       projection angles 
%   ++par         parameter structure 
% 

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
%   using the “cSAXS tomography package” developed by the CXS group,
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


function [stack_object, theta, par] = load_stored_object(path)
    utils.verbose(0, 'Loading cached object ..... ')
    load(path)
    try
        % create a complex number 
        stack_object = complex(stack_object_r, stack_object_i); 
    catch
        try
            stack_object_r = complex(stack_object_r);
            % create a complex number inplace of the stack_object_r
            stack_object = tomo.block_fun(@(re,im)(re+1i*im),stack_object_r, stack_object_i, struct('use_GPU',false, 'inplace', true)); 
        catch err
            warning('Loading failed')
            disp(err)
            keyboard
        end
    end
    utils.verbose(0, 'Loading cached object done ')
end