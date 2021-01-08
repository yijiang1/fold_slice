/*

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
   proper use and the correctness of the results.


*/


// Defines the exported functions for the DLL application.
//
// recompile commands
//  (Linux, GCC 4.8.5)   mexcuda -outdir private  ASTRA_GPU_wrapper/ASTRA_GPU_wrapper.cu ASTRA_GPU_wrapper/util3d.cu ASTRA_GPU_wrapper/par3d_fp.cu ASTRA_GPU_wrapper/par3d_bp.cu
//  (Windows)  mexcuda -outdir private  ASTRA_GPU_wrapper\ASTRA_GPU_wrapper.cu ASTRA_GPU_wrapper\util3d.cu ASTRA_GPU_wrapper\par3d_fp.cu ASTRA_GPU_wrapper\par3d_bp.cu

/************* INPUTS *****************************/
/* 
    string 'fp' or 'bp'  - forward / backward projection 
    single gpuArray volume or data object 
    struct cfg - contain configuration for astra, created by ASTRA_initialize.m
    double array vec - contain projection geometry for astra, created by ASTRA_initialize.m
    (optional)
    single gpuArray - volume or data object to write the results to 
*/


#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda.h>

#include <stdio.h>

#include <cstdio>
#include <cassert>
#include <iostream>
#include <list>

#include "mex.h"
#include "gpu/mxGPUArray.h"

#include "util3d.h"
#include "dims3d.h"
#include "par3d_bp.h"
#include "par3d_fp.h"



void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, mxArray const *prhs[])
{

    //mexPrintf("Warning: loading development version of ASTRA\n");
    //mexPrintf("Ninputs:%i\n", nrhs); 

    if (!((nrhs == 4) || (nrhs == 5) || (nrhs == 8 ) || (nrhs == 11 ) ))
        mexErrMsgTxt("4,5, 8, or 11 input arguments required");


	using namespace astraCUDA3d;
	char const * const errId = "parallel:gpu:mexGPUExample:InvalidInput";
	char const * const errMsg = "Invalid input to MEX file.";


	/* Throw an error if the input is not a GPU array. */
	if (!mxIsGPUArray(prhs[1])) {
		mexErrMsgIdAndTxt(errId, "The second  input must be GPU array");
	}
	

	/* Load configuration */
	SDimensions3D dims;
	mxArray * tmp; 
	double * val;
    #define SETVAR(name) do {tmp = mxGetField(prhs[2], 0, ""#name"");  if (tmp!=NULL) { val = mxGetPr(tmp); dims.name = (unsigned int)val[0]; }} 	while (0); 
        SETVAR(iVolX);
        SETVAR(iVolY);
        SETVAR(iVolZ);
        SETVAR(iProjAngles);
        SETVAR(iProjU);
        SETVAR(iProjV);
        SETVAR(iRaysPerDetDim);
        SETVAR(iRaysPerVoxelDim);
    #undef SETVAR


	/* Initialize the MathWorks GPU API. */
	mxInitGPU();


	/* load confuguration of angles */
	double * my_angles = mxGetPr(prhs[3]);
	int Nangles = (int)mxGetM(prhs[3]);
	SPar3DProjection* angle = new SPar3DProjection[Nangles];

#define SETVAR(name,i,j) do { angle[i].name = my_angles[i+j*Nangles];  } while (0); 
	for (int i = 0; i < Nangles; i++)
	{
		SETVAR(fRayX, i, 0);
		SETVAR(fRayY, i, 1);
		SETVAR(fRayZ, i, 2);
		SETVAR(fDetSX, i, 3);
		SETVAR(fDetSY, i, 4);
		SETVAR(fDetSZ, i, 5);
		SETVAR(fDetUX, i, 6);
		SETVAR(fDetUY, i, 7);
		SETVAR(fDetUZ, i, 8);
		SETVAR(fDetVX, i, 9);
		SETVAR(fDetVY, i, 10);
		SETVAR(fDetVZ, i, 11);
		// mexPrintf("---------------------- \n"); 
	}
#undef SETVAR


	char * task = mxArrayToString(prhs[0]);

	//mexPrintf("--------- Task %s \n ", task); 


	/* Load input data */
	mxGPUArray const * m_data = mxGPUCreateFromMxArray(prhs[1]);
	if ((mxGPUGetClassID(m_data) != mxSINGLE_CLASS)) {
		mexErrMsgIdAndTxt(errId, errMsg);
	}
	float * p_data = (float *)mxGPUGetDataReadOnly(m_data);

    DeformField DF;
    if (nrhs == 8 || nrhs == 11 ) {
        /* load deformation field */
        DF.use_deform = true;
        DF.use_linear_model = false; // assume contant deformation 
        DF.X0 = mxGPUCreateFromMxArray(prhs[5]);
        DF.Y0 = mxGPUCreateFromMxArray(prhs[6]);
        DF.Z0 = mxGPUCreateFromMxArray(prhs[7]);
        if ((mxGPUGetClassID(DF.X0) != mxSINGLE_CLASS) |
            (mxGPUGetClassID(DF.Y0) != mxSINGLE_CLASS) |
            (mxGPUGetClassID(DF.Z0) != mxSINGLE_CLASS)) {
           mexPrintf("wrong input type: deformation fields has to be single\n");
           mexErrMsgIdAndTxt(errId, errMsg);
        }
        if (nrhs == 11 ) {
            DF.use_linear_model = true; // assume linear deformation 
            DF.X1 = mxGPUCreateFromMxArray(prhs[8]);
            DF.Y1 = mxGPUCreateFromMxArray(prhs[9]);
            DF.Z1 = mxGPUCreateFromMxArray(prhs[10]);
            if ((mxGPUGetClassID(DF.X1) != mxSINGLE_CLASS) |
                (mxGPUGetClassID(DF.Y1) != mxSINGLE_CLASS) |
                (mxGPUGetClassID(DF.Z1) != mxSINGLE_CLASS)) {
               mexPrintf("wrong input type: deformation fields has to be single\n");
               mexErrMsgIdAndTxt(errId, errMsg);
            }
        }
        }
    else 
        DF.use_deform = false;
  


	if (strcmp(task, "fp")==0)
	{
		//mexPrintf(" forward projection \n ");

		/* make volume array (no copying) */
		cudaPitchedPtr volData;
		volData.ptr = p_data; 
		volData.pitch = dims.iVolX * sizeof(float);
		volData.xsize = dims.iVolX;
		volData.ysize = dims.iVolY;


    mxGPUArray * m_projData;
    if(nrhs >= 5 && !mxIsEmpty(prhs[4]) ) 
    {
        /**** copy of the array is the slow operation and also GPU memory is limited *****/ 
       // m_projData = mxGPUCopyFromMxArray(prhs[4]); 

        /*  Use ugly trick to write directly to the provided GPU array ... 
          =>  Now it is writting directly into the input field !!! DANGEROUS  */

        m_projData = const_cast<mxGPUArray*>(mxGPUCreateFromMxArray(prhs[4])); 
        if ((mxGPUGetClassID(m_projData) != mxSINGLE_CLASS)) {
               mexPrintf("m_projData\n");
               mexErrMsgIdAndTxt(errId, errMsg);
         }
        const mwSize * projSize = mxGPUGetDimensions(m_projData);

        if (dims.iProjU != projSize[0]   || 
            dims.iProjV != projSize[1]    ||
            dims.iProjAngles != projSize[2])
            mexErrMsgIdAndTxt(errId, "Wrong size of the inputs array");

        //mexPrintf("Writting directly to the input array\n\n");

    }
    else
    {
		/* allocate projection field  */
		int const Ndim = 3;
		mwSize projSize[3];
		projSize[0] = (mwSize)dims.iProjU;
		projSize[1] = (mwSize)dims.iProjV;
		projSize[2] = (mwSize)dims.iProjAngles;
		m_projData = mxGPUCreateGPUArray(Ndim,
			projSize,
			mxSINGLE_CLASS,
			mxREAL,
			MX_GPU_INITIALIZE_VALUES);
    }

		/* make cudaPitchedPtr for  projection field  */
		cudaPitchedPtr projData;
		projData.ptr = (float *)mxGPUGetData(m_projData);
		projData.pitch = dims.iProjU * sizeof(float);
		projData.xsize = dims.iProjU;
		projData.ysize = dims.iProjV;

        //mexPrintf("astraCUDA3d::Par3DFP \n ") ;
		astraCUDA3d::Par3DFP(volData, projData, dims, angle, 1.0f, DF);
		checkLastError("After Projector");


		/* Wrap the result up as a MATLAB gpuArray for return. */
        if (nlhs > 0)
    		plhs[0] = mxGPUCreateMxArrayOnGPU(m_projData);
		mxGPUDestroyGPUArray(m_projData);
		mxGPUDestroyGPUArray(m_data);


	}
	else if (strcmp(task, "bp")==0)
	{
		//mexPrintf(" backward projection \n  ");


		/* make  projection field (no copying) */
		cudaPitchedPtr projData;
		projData.ptr = p_data;
		projData.pitch = dims.iProjU * sizeof(float);
		projData.xsize = dims.iProjU;
		projData.ysize = dims.iProjAngles;
        mxGPUArray* m_volData; 
        if(nrhs >= 5 && !mxIsEmpty(prhs[4]) ) 
        {
            /**** copy of the array is the slow operation and also GPU memory is limited *****/ 
           // m_volData = mxGPUCopyFromMxArray(prhs[4]); 

            /*  Use ugly trick to write directly to the provided GPU array ... 
              =>  Now it is writting directly into the input field !!! DANGEROUS  */

            m_volData = const_cast<mxGPUArray*>(mxGPUCreateFromMxArray(prhs[4])); 
            if ((mxGPUGetClassID(m_volData) != mxSINGLE_CLASS)) {
                   mexPrintf("m_volData\n");
                   mexErrMsgIdAndTxt(errId, errMsg);
             }
            mwSize volSize[3];
            const mwSize * volSize0 = mxGPUGetDimensions(m_volData);
            if (mxGPUGetNumberOfDimensions(m_volData)==3) {
                volSize[0]=volSize0[0];
                volSize[1]=volSize0[1];
                volSize[2]=volSize0[2];
            } else {
                volSize[0]=volSize0[0];
                volSize[1]=volSize0[1];
                volSize[2]=1;
            }


            if (dims.iVolX != volSize[0]   || 
                dims.iVolY != volSize[1]    ||
                dims.iVolZ != volSize[2])
                mexErrMsgIdAndTxt(errId, "Wrong size of the inputs array");

        } else {
            /* allocate volume data  */
            int const Ndim = 3;
            mwSize volSize[3];
            volSize[0] = (mwSize)dims.iVolX;
            volSize[1] = (mwSize)dims.iVolY;
            volSize[2] = (mwSize)dims.iVolZ;
            m_volData = mxGPUCreateGPUArray(Ndim,
                volSize,
                mxSINGLE_CLASS,
                mxREAL,
                MX_GPU_INITIALIZE_VALUES);
        }

		/* make volume array pointer*/
		cudaPitchedPtr volData;
		volData.ptr = (float *)mxGPUGetData(m_volData);
		volData.pitch = dims.iVolX * sizeof(float);
		volData.xsize = dims.iVolX;
		volData.ysize = dims.iVolY;

		astraCUDA3d::Par3DBP(volData, projData, dims, angle, 1.0f, DF);
		checkLastError("After Projector");


		/* Wrap the result up as a MATLAB gpuArray for return. */
        if (nlhs > 0)
            plhs[0] = mxGPUCreateMxArrayOnGPU(m_volData);
		mxGPUDestroyGPUArray(m_volData);
		mxGPUDestroyGPUArray(m_data);

	}
	else
		mexPrintf("No such option"); 

    if (DF.use_deform) {
        //mexPrintf("Deleted DF"); 
        mxGPUDestroyGPUArray(DF.X0);
        mxGPUDestroyGPUArray(DF.Y0);
        mxGPUDestroyGPUArray(DF.Z0);
        if (DF.use_linear_model) {
            mxGPUDestroyGPUArray(DF.X1);
            mxGPUDestroyGPUArray(DF.Y1);
            mxGPUDestroyGPUArray(DF.Z1);
        }
    }


}

