#include "mex.h"
#include "gpu/mxGPUArray.h"
#include "interp3_gpu.hpp"

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


/**
 * MEX gateway
 */
void mexFunction(int nlhs , mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{
    char const * const errId = "parallel:gpu:interp3_gpu:InvalidInput";
    char const * const errMsg = "Invalid input to MEX file.";

    // Initialize the MathWorks GPU API.
    mxInitGPU();

    if (nrhs!=4) {
        mexPrintf("Wrong number of inputs\n");
        mexErrMsgIdAndTxt(errId, errMsg);
    }
    
   
    const mxGPUArray * m_Img_orig = mxGPUCreateFromMxArray(prhs[0]); 
	if ((mxGPUGetClassID(m_Img_orig) != mxSINGLE_CLASS)) {
       mexPrintf("wrong input m_Img_orig\n");
	   mexErrMsgIdAndTxt(errId, errMsg);
	}
	const float * p_Img_orig = (const float *)mxGPUGetDataReadOnly(m_Img_orig);
 
    
    const mxGPUArray * m_X = mxGPUCreateFromMxArray(prhs[1]);
    const mxGPUArray * m_Y = mxGPUCreateFromMxArray(prhs[2]);
    const mxGPUArray * m_Z = mxGPUCreateFromMxArray(prhs[3]);
    if ((mxGPUGetClassID(m_X) != mxSINGLE_CLASS) |
        (mxGPUGetClassID(m_Y) != mxSINGLE_CLASS) |
        (mxGPUGetClassID(m_Z) != mxSINGLE_CLASS)) {
       mexPrintf("wrong input X,Y,Z\n");
	   mexErrMsgIdAndTxt(errId, errMsg);
	}

    mwSize  const * dimensions = mxGPUGetDimensions(m_Img_orig);
    mwSize Ndim =  mxGPUGetNumberOfDimensions(m_Img_orig);
    int M = (int)dimensions[0];
    int N = (int)dimensions[1];
    int O = Ndim > 2 ? (int)dimensions[2] :  1;
            
    mxGPUArray * m_Img_out = mxGPUCreateGPUArray(Ndim,
                dimensions,
                mxSINGLE_CLASS,
                mxREAL,
                MX_GPU_INITIALIZE_VALUES);
    float *  p_Img_out = (float *)mxGPUGetData(m_Img_out);
    
    checkLastError("Before kernel run");

    // mexcuda -output interp3_gpu_mex interp3_gpu_ker.cu interp3_gpu_mex.cpp
        
    interp3_init( p_Img_out,m_Img_orig, m_X, m_Y, m_Z, M, N, O);
 
   checkLastError("Before after run");

    plhs[0] = mxGPUCreateMxArrayOnGPU(m_Img_out);

    mxGPUDestroyGPUArray(m_Img_out);
    mxGPUDestroyGPUArray(m_Img_orig);
    mxGPUDestroyGPUArray(m_X);
    mxGPUDestroyGPUArray(m_Y);
    mxGPUDestroyGPUArray(m_Z);
    
    
}


