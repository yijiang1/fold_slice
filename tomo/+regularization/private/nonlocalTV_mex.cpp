
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include "TV_texture.hpp"

/**
 *
 **-----------------------------------------------------------------------*
|                                                                       |
|  Except where otherwise noted, this work is licensed under a          |
|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
|  International (CC BY-NC-SA 4.0) license.                             |
|                                                                       |
|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
|                                                                       |
|       Author: CXS group, PSI                                          |
*-----------------------------------------------------------------------*
You may use this code with the following provisions:

If the code is fully or partially redistributed, or rewritten in another
  computing language this notice should be included in the redistribution.

If this code, or subfunctions or parts of it, is used for research in a 
  publication or if it is fully or partially rewritten for another 
  computing language the authors and institution should be acknowledged 
  in written form in the publication: “Data processing was carried out 
  using the “cSAXS matlab package” developed by the CXS group,
  Paul Scherrer Institut, Switzerland.” 
  Variations on the latter text can be incorporated upon discussion with 
  the CXS group if needed to more specifically reflect the use of the package 
  for the published work.

A publication that focuses on describing features, or parameters, that
   are already existing in the code should be first discussed with the
   authors.
  
This code and subroutines are part of a continuous development, they 
   are provided “as they are” without guarantees or liability on part
   of PSI or the authors. It is the user responsibility to ensure its 

 *
 *
 * MEX gateway
 */
void mexFunction(int nlhs , mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{
    char const * const errId = "parallel:gpu:mexGPUExample:InvalidInput";
    char const * const errMsg = "Invalid input to MEX file.";

    // Initialize the MathWorks GPU API.
    mxInitGPU();

    if (nrhs!=10) {
       mexPrintf("nargin\n");
        mexErrMsgIdAndTxt(errId, errMsg);
    }
    
    // Im0,CloseInd,neighbours, dt, eps, lambda, Nclose, Rwin, Niter, Nx, Ny

        
        
    mxGPUArray const * m_Img = mxGPUCreateFromMxArray(prhs[0]);
	if ((mxGPUGetClassID(m_Img) != mxSINGLE_CLASS)) {
       mexPrintf("m_Img\n");
	   mexErrMsgIdAndTxt(errId, errMsg);
	}
	const float * p_Img = (float *)mxGPUGetDataReadOnly(m_Img);
   
    mxGPUArray const * m_CloseInd = mxGPUCreateFromMxArray(prhs[1]);
	if ((mxGPUGetClassID(m_CloseInd) != mxUINT8_CLASS)) {
       mexPrintf("m_CloseInd\n");
	   mexErrMsgIdAndTxt(errId, errMsg);
	}
	const uint8_T * p_CloseInd = (uint8_T *)mxGPUGetDataReadOnly(m_CloseInd);
    

    mxGPUArray const * m_Nbrs = mxGPUCreateFromMxArray(prhs[2]);
	if ((mxGPUGetClassID(m_Nbrs) != mxUINT8_CLASS)) {
       mexPrintf("m_Nbrs\n");
	   mexErrMsgIdAndTxt(errId, errMsg);
	}
	const uint8_T * p_Nbrs = (uint8_T *)mxGPUGetDataReadOnly(m_Nbrs);

    mxGPUArray const * mNbrs_weights = mxGPUCreateFromMxArray(prhs[3]);
    if ((mxGPUGetClassID(mNbrs_weights) != mxSINGLE_CLASS)) {
       mexPrintf("mNbrs_weights\n");
	   mexErrMsgIdAndTxt(errId, errMsg);
	}
	const float * p_Nbrs_weights = (float *)mxGPUGetDataReadOnly(mNbrs_weights);

    
    
   const float dt = (float)mxGetScalar(prhs[4]);
   const float eps = (float)mxGetScalar(prhs[5]);
   const float lambda = (float)mxGetScalar(prhs[6]);

   const int Nclose = (int)mxGetScalar(prhs[7]);
   const int Rwin = (int)mxGetScalar(prhs[8]);
   
   const int Niter = (int)mxGetScalar(prhs[9]);
   
   const int Nnbrs = mxGetNumberOfElements(prhs[2]); 
   

    mwSize  const * dimensions = mxGPUGetDimensions(m_Img);
    mwSize Ndim =  mxGPUGetNumberOfDimensions(m_Img);
    int M = (int)dimensions[0];
    int N = (int)dimensions[1];
    int O = Ndim > 2 ? (int)dimensions[2] :  1;
    

 
    /* allocate output image field  */
    mxGPUArray * m_Img_new = mxGPUCopyFromMxArray(prhs[0]); 
	float * p_Img_new = (float *)mxGPUGetData(m_Img_new);
    
    
//     for(int i = 0;  i < Nnbrs; i++)
//         mexPrintf("neighbours %i\n", p_Nbrs[i]);
//    for(int i = 0;  i < Nnbrs; i++)
//        mexPrintf("Nbrs_weights %g\n", p_Nbrs_weights[i]);


   
    nonlocal_TV_init( p_Img_new, p_Img, p_CloseInd, p_Nbrs, p_Nbrs_weights, dt, eps,
                        lambda, Nclose, Nnbrs,  Rwin, N,  M, O, Niter);
    checkLastError("After iteration");
   
 
   
    plhs[0] = mxGPUCreateMxArrayOnGPU(m_Img_new);
    mxGPUDestroyGPUArray(m_CloseInd);
    mxGPUDestroyGPUArray(m_Img);   
    mxGPUDestroyGPUArray(m_Nbrs);
    mxGPUDestroyGPUArray(m_Img_new);
    mxGPUDestroyGPUArray(mNbrs_weights);

    
    
}
