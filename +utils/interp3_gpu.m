% INTERP3_GPU - fast texture-based GPU based interpolation method for 3D deformation 
% input array is deformated gived X,Y,Z deformation vector fields 
%
% array = interp3_gpu(array, DVF_X, DVF_Y, DVF_Z)
%
% Inputs:
%    **array            volume to be deformed 
%    **DVF_X            deformation field in X direction 
%    **DVF_Y            deformation field in Y direction 
%    **DVF_Z            deformation field in Z direction 
% Outputs: 
%    ++array            deformed volume 


%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2018 by Paul Scherrer Institute (http://www.psi.ch)    |
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


function array = interp3_gpu(array, DVF_X, DVF_Y, DVF_Z)


  % apply 3D deformation using GPU textures 
  try
     array = interp3_gpu(array, DVF_X, DVF_Y, DVF_Z); 
  catch err 
      if strcmpi(err.identifier, 'MATLAB:mex:ErrInvalidMEXFile')
            % recompile the MEX code 
            path = replace(mfilename('fullpath'), mfilename, ''); 
            mexcuda('-output', fullfile(path,'private/interp3_gpu_ker'), fullfile(path, 'private/interp3_gpu_ker.cu'))
            array = interp3_gpu(array, DVF_X, DVF_Y, DVF_Z);    
      else
          rethrow(err) 
      end
  end
      


end