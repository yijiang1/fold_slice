
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

    if (nrhs!=7) {
       mexPrintf("nargin\n");
        mexErrMsgIdAndTxt(errId, errMsg);
    }
    
    mxGPUArray const * m_Img = mxGPUCreateFromMxArray(prhs[0]);
	if ((mxGPUGetClassID(m_Img) != mxSINGLE_CLASS)) {
       mexPrintf("m_Img\n");
       mexErrMsgIdAndTxt(errId, errMsg);
	}
	const float * p_Img = (float *)mxGPUGetDataReadOnly(m_Img);
    
    

    mxGPUArray const * m_Nbrs = mxGPUCreateFromMxArray(prhs[1]);
	if ((mxGPUGetClassID(m_Nbrs) != mxUINT8_CLASS)) {
       mexPrintf("m_Nbrs\n");

		mexErrMsgIdAndTxt(errId, errMsg);
	}
	const uint8_T * p_Nbrs = (uint8_T *)mxGPUGetDataReadOnly(m_Nbrs);

    
   const int Nclose = (int)mxGetScalar(prhs[2]);
   const int Nclose_min = (int)mxGetScalar(prhs[3]);

   const int Rwin = (int)mxGetScalar(prhs[4]);
   const int Rpatch = (int)mxGetScalar(prhs[5]);
   const float threshold = (float)mxGetScalar(prhs[6]);

   //mexPrintf("thresh %3.5g \n ", threshold); 

   

    mwSize  const * dimensions = mxGPUGetDimensions(m_Img);
    mwSize const Ndim =  mxGPUGetNumberOfDimensions(m_Img);
    const int M = (int)dimensions[0];
    const int N = (int)dimensions[1];
    const int O = Ndim > 2 ? (int)dimensions[2] :  1;
    
   // mexPrintf("%i %i %i \n ", M,N,O); 
   // mexPrintf("Nc %i Rw %i Rp%i \n ", Nclose, Rwin, Rpatch); 

   const int Nnbrs = mxGetNumberOfElements(prhs[1]); 


  
    /* allocate index of closest patches field  */
    mwSize matSize[4];
    matSize[0] = M;
    matSize[1] = N;
    matSize[2] = O;
    matSize[3] = Nclose;
    mxGPUArray* m_CloseInd = mxGPUCreateGPUArray(4,
        matSize,
        mxUINT8_CLASS,
        mxREAL,
        MX_GPU_INITIALIZE_VALUES);
    
    
    uint8_T * p_CloseInd = (uint8_T *)mxGPUGetData(m_CloseInd);

  //  mexPrintf("%i %i %i \n ", matSize[0], matSize[1], matSize[2]); 

    

    nonlocal_weight_TV_init( p_CloseInd, p_Img, p_Nbrs,  Nclose,Nclose_min, Nnbrs,  Rwin, Rpatch, N,  M, O, threshold);
   
 
   
    plhs[0] = mxGPUCreateMxArrayOnGPU(m_CloseInd);
    mxGPUDestroyGPUArray(m_CloseInd);
    mxGPUDestroyGPUArray(m_Img);   
    mxGPUDestroyGPUArray(m_Nbrs);

    
    
}
