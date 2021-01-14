#ifndef INTERP3_GPU_tex_HPP
#define INTERP3_GPU_tex_HPP

#include "tmwtypes.h"
#include "mex.h"
#include "gpu/mxGPUArray.h"

//  interp3_gpu.m - fast texture based GPU based interpolation method for 3D deformation 
//
//  RECOMPILE: mexcuda -output interp3_gpu_mex interp3_gpu_ker.cu interp3_gpu_mex.cpp
// 
// %*-----------------------------------------------------------------------*
// %|                                                                       |
// %|  Except where otherwise noted, this work is licensed under a          |
// %|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
// %|  International (CC BY-NC-SA 4.0) license.                             |
// %|                                                                       |
// %|  Copyright (c) 2018 by Paul Scherrer Institute (http://www.psi.ch)    |
// %|                                                                       |
// %|       Author: CXS group, PSI                                          |
// %*-----------------------------------------------------------------------*
// % You may use this code with the following provisions:
// %
// % If the code is fully or partially redistributed, or rewritten in another
// %   computing language this notice should be included in the redistribution.
// %
// % If this code, or subfunctions or parts of it, is used for research in a 
// %   publication or if it is fully or partially rewritten for another 
// %   computing language the authors and institution should be acknowledged 
// %   in written form in the publication: “Data processing was carried out 
// %   using the “cSAXS matlab package” developed by the CXS group,
// %   Paul Scherrer Institut, Switzerland.” 
// %   Variations on the latter text can be incorporated upon discussion with 
// %   the CXS group if needed to more specifically reflect the use of the package 
// %   for the published work.
// %
// % A publication that focuses on describing features, or parameters, that
// %    are already existing in the code should be first discussed with the
// %    authors.
// %   
// % This code and subroutines are part of a continuous development, they 
// %    are provided “as they are” without guarantees or liability on part
// %    of PSI or the authors. It is the user responsibility to ensure its 
// %    proper use and the correctness of the results.




int checkLastError(char * msg);

void interp3_init( float * Img, const mxGPUArray * Img_0, const mxGPUArray *X, const mxGPUArray *Y, const mxGPUArray *Z, const unsigned int M, const unsigned int N, const unsigned int O);


#endif

