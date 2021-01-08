#include "mex.h"
#include "gpu/mxGPUArray.h"
#include "TV_texture.hpp"

    //mexcuda -output private/local_TV_mex private/TV_cuda_texture.cu private/local_TV_mex.cpp


/**
 *


*-----------------------------------------------------------------------*
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
 * MEX gateway
 */
void mexFunction(int nlhs , mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{
    char const * const errId = "parallel:gpu:mexGPUExample:InvalidInput";
    char const * const errMsg = "Invalid input to MEX file.";

    // Initialize the MathWorks GPU API.
    mxInitGPU();

    if (nrhs!=5) {
        mexPrintf("Wrong number of inputs\n");
        mexErrMsgIdAndTxt(errId, errMsg);
    }
    
   
    
   const float dt = (float)mxGetScalar(prhs[1]);
   const float tau = (float)mxGetScalar(prhs[2]);
   const int Niter = (int)mxGetScalar(prhs[3]);
   const int UseChambolle = (int)mxGetScalar(prhs[4]);

    /* allocate output image field  */
    mxGPUArray * m_Img_new = mxGPUCopyFromMxArray(prhs[0]); 
	if ((mxGPUGetClassID(m_Img_new) != mxSINGLE_CLASS)) {
       mexPrintf("m_Img_new\n");
	   mexErrMsgIdAndTxt(errId, errMsg);
	}
	float * p_Img_new = (float *)mxGPUGetData(m_Img_new);
 
 
    mwSize  const * dimensions = mxGPUGetDimensions(m_Img_new);
    mwSize Ndim =  mxGPUGetNumberOfDimensions(m_Img_new);
    int M = (int)dimensions[0];
    int N = (int)dimensions[1];
    int O = Ndim > 2 ? (int)dimensions[2] :  1;
    mxGPUArray * mXi[3] ;

    if (UseChambolle ) {
        /* allocate Xi field  */
        float * Xi[3] ;
        mwSize xiSize[3] = {M,N,O};

        for( int i = 0; i < 3; i++) {
            mXi[i]  = mxGPUCreateGPUArray(3,
                xiSize,
                mxSINGLE_CLASS,
                mxREAL,
                MX_GPU_INITIALIZE_VALUES);
           Xi[i] = (float *)mxGPUGetData(mXi[i]);
        }
        //mexPrintf("local_TV_chambolle_init\n");
        local_TV_chambolle_init( p_Img_new, Xi, dt, tau, M, N, O, Niter);
    } else {
       // mexPrintf("local_TV_init\n");
         local_TV_init( p_Img_new, dt, tau,  N,  M, O, Niter);       
    }
    
    
            
    checkLastError("After iteration");
   
    //mexcuda -output utils3D\local_TV_mex utils3D\TV_cuda_texture.cu utils3D\local_TV_mex.cpp
            
 
   
    plhs[0] = mxGPUCreateMxArrayOnGPU(m_Img_new);

    mxGPUDestroyGPUArray(m_Img_new);
    if (UseChambolle )
        for( int i = 0; i < 3; i++)
            mxGPUDestroyGPUArray( mXi[i]);

    
    
}

