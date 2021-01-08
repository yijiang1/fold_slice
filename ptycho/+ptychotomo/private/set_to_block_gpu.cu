/*  
    SET_TO_BLOCK_GPU  Set complex views to a large complex object without memory copy of the large volume 

    set_to_block_gpu(volData, volData_block, uint16(ind));        

  Inputs: 
   **volData         - 3D array, single precision 
   **volData_block   - 3D array, that will be written into the volData
   **ind             - 1D array, uint16, indices into which lines of volData will be volData_block written 
  Outputs:
    none, results are written to volData directly 

 Recompile: 
    mexcuda -output +engines/+GPU/+shared/private/set_views_gpu_mex +engines/+GPU/+shared/private/set_views_gpu_mex.cu
 */

#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <math.h> 
#include <stdio.h> 

typedef const unsigned int   cuint;
typedef const uint16_T cuint16;

// unfortunatelly ~10800 is the maximum of const memory 
const unsigned int MAX_IND_READ = 10800;
static const unsigned MAX_IND_READ_DEV = MAX_IND_READ;
__constant__ uint16_T gC_ind[MAX_IND_READ_DEV];


int checkLastError(char * msg)
{
    cudaError_t cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        char err[512];
        sprintf(err, "setprojection failed \n %s: %s. \n", msg, cudaGetErrorString(cudaStatus));
        mexPrintf(err);
        return 1;
    }
    return 0;
}


/*
 * Device code
 */


__global__ void addToArray_c( float2 const * sarray, float2 * larray, cuint Nx,cuint Ny, cuint Npos)  {
    // Location in a 3D matrix
    int idx= blockIdx.x * blockDim.x + threadIdx.x;
    int idy= blockIdx.y * blockDim.y + threadIdx.y;
    if ( idx < Nx & idy < Ny ) {
        int idz, id_large, id_small; 
        for(idz = 0; idz < Npos; idz++)
        {   
            id_large = idx + Nx*idy  + Nx*Ny*(gC_ind[idz]-1);  // go only through some of the indices 
            id_small = idx + Nx*idy  + Nx*Ny*idz; 


            larray[id_large].x  = sarray[ id_small ].x;
            larray[id_large].y  = sarray[ id_small ].y;
            
        }
    }
}


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{


    int i ; 
    char const * const errId = "parallel:gpu:mexGPUExample:InvalidInput";
    char const * const errMsg = "Invalid input to MEX file.";
    
    // Check for proper number of arguments.
    if (nrhs != 3) 
    mexErrMsgTxt("Three input arguments required");
  

    // Input must be of type single. 
  for (i=0; i < 2; i++) {
      if ( !mxIsGPUArray(prhs[i]) ){
          printf("Input %d is not cell array\n",i+1);
           mexErrMsgIdAndTxt("MexError:ptycho","Inputs must be of correct type.");
      }
  }
    // Input must be of type int16.
  for (i=2; i<3; i++){
    if (mxIsUint16(prhs[i]) != 1){
          printf("Input %d is not integer\n",i+1);
          mexErrMsgIdAndTxt("MexError:ptycho","Inputs must be of correct type uint16.");
    }
  }
 
   // It cannot be one-dimensional 
  if(mxGetNumberOfDimensions(prhs[0]) > 3) {
    printf("The 1st input argument must have at least three dimensions.");
    mexErrMsgIdAndTxt("MexError:ptycho","wrong number of dimensions");
  }
    // It cannot be more than 3-dimensional 
   if(mxGetNumberOfDimensions(prhs[0]) > 3) {
    printf("The 1st input argument must have at most three dimensions.");
    mexErrMsgIdAndTxt("MexError:ptycho","wrong number of dimensions");
  }

  // Check that arrays are complex


     const mxGPUArray * m_sarray = mxGPUCreateFromMxArray(prhs[1]); 
     if ((mxGPUGetClassID(m_sarray) != mxSINGLE_CLASS)) {
           mexPrintf("m_sarray\n");
           mexErrMsgIdAndTxt(errId, errMsg);
     }

    // write results directly to the output array 
     mxGPUArray * m_larray = const_cast<mxGPUArray*>(mxGPUCreateFromMxArray(prhs[0])); //    mxGPUCopyFromMxArray(prhs[1]); 
     if ((mxGPUGetClassID(m_larray) != mxSINGLE_CLASS)) {
           mexPrintf("m_larray\n");
           mexErrMsgIdAndTxt(errId, errMsg);
     }

     if (mxGPUGetComplexity(m_larray) != mxCOMPLEX) {
        mexPrintf("m_larray complexity has to be the COMPLEX \n");
           mexErrMsgIdAndTxt(errId, errMsg);
     }


     if (mxGPUGetComplexity(m_larray) != mxGPUGetComplexity(m_sarray)) {
        mexPrintf("m_larray and m_sarray complexity has to be the same \n");
           mexErrMsgIdAndTxt(errId, errMsg);
     }



     const mxGPUArray * m_ind = mxGPUCreateFromMxArray(prhs[2]); 
     if ((mxGPUGetClassID(m_ind) != mxUINT16_CLASS)) {
           mexPrintf("m_ind\n");
           mexErrMsgIdAndTxt(errId, errMsg);
     }
     const uint16_T *  p_ind = (uint16_T *)mxGPUGetDataReadOnly(m_ind);



    // Get dimension of probe and object 
    const unsigned int  Ndims = (unsigned int)mxGPUGetNumberOfDimensions(m_sarray);
    const mwSize * Np_l = mxGPUGetDimensions(m_larray);
    const mwSize * Np_s = mxGPUGetDimensions(m_sarray);
    const unsigned int Npos = mxGPUGetNumberOfElements(m_ind);




    if (Ndims == 3 && Npos > Np_s[2]) {
      printf("wrong size of update / positions  %i", Ndims);
      mexErrMsgIdAndTxt("MexError:ptycho","wrong size of update / positions");
    } 


    cudaMemcpyToSymbol(gC_ind, p_ind, Npos*sizeof(uint16_T), 0, cudaMemcpyHostToDevice);
    checkLastError("after cudaMemcpyToSymbol pos");


    // Choose a reasonably sized number of threads in each dimension for the block.
    int const threadsPerBlockEachDim =  32;
    // Compute the thread block and grid sizes based on the board dimensions.
    int const blocksPerGrid_M = (Np_s[0] + threadsPerBlockEachDim - 1) / threadsPerBlockEachDim;
    int const blocksPerGrid_N = (Np_s[1] + threadsPerBlockEachDim - 1) / threadsPerBlockEachDim;
    int const blocksPerGrid_O = 1;

    dim3 const dimBlock(blocksPerGrid_M, blocksPerGrid_N, blocksPerGrid_O);
    dim3 const dimThread(threadsPerBlockEachDim, threadsPerBlockEachDim, 1);

    checkLastError("after dimThread");

    checkLastError("after cudaMemcpyToSymbol");

    // ================== call the right kernel =================== 

    const float2 * p_sarray = (float2 *)mxGPUGetDataReadOnly(m_sarray);
    float2 * p_larray = (float2 *)mxGPUGetData(m_larray);
    addToArray_c<<<dimBlock, dimThread>>>(p_sarray,p_larray, Np_l[0],Np_l[1], Npos);


    checkLastError("after kernel");

    mxGPUDestroyGPUArray(m_larray);
    mxGPUDestroyGPUArray(m_ind);

    cudaThreadSynchronize();

    // plhs[0] = mxGPUCreateMxArrayOnGPU(m_larray);
    mxGPUDestroyGPUArray(m_sarray);

    
  return;
}




