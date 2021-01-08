% SAVE_TOMOGRAM simple function to save data as tiff images 
% save current tomogram 
%
% save_tomogram(tomogram, par, type, circulo,theta, extra_string = '')
%
% Inputs:
%   **tomogram          - (3D array) saved volume 
%   **par               - tomo parameter structure 
%   **type              - 'delta', 'beta'
%   **circulo, theta    - other inputs to be saved 
%   **extra_string      - extra string in the name, default == ''

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


function save_tomogram(tomogram, par, type, circulo,theta, extra_string)
    % save current tomogram 
    % type: 'delta', 'beta', ''
    
    if nargin < 6
        extra_string = ''; 
    end
    if ~exist(par.output_folder,'dir')
        mkdir(par.output_folder); 
    end
      

    saveprojfile = fullfile(sprintf('%s/tomogram_%s_%s_%s.mat',par.output_folder,type,par.scans_string, extra_string));
    
    switch type
        case 'delta'
            tomogram_delta = tomogram; 
        case 'beta'
            tomogram_beta = tomogram; 
        otherwise
            error('Allowed types: "delta", "beta"')
    end
    
    utils.verbose(0, 'Saving tomogram ....')
    utils.savefast_safe(saveprojfile,['tomogram_',type],'par','circulo', 'theta', par.force_overwrite);
    utils.verbose(0, 'Saving done')
    
end
