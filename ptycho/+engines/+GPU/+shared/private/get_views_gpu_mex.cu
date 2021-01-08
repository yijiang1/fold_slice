/*  
    Get complex views from complex object 
    mexcuda -output +engines/+GPU/+shared/private/get_views_gpu_mex +engines/+GPU/+shared/private/get_views_gpu_mex.cu
 */

#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <math.h> 
#include <stdio.h> 
typedef const unsigned int   cuint;
typedef const uint16_T cuint16;

// unfortunatelly ~10800 is the maximum of const memory 
const unsigned int MAX_IND_READ = 10800;
__constant__ uint16_T gC_ind_read[MAX_IND_READ];
__constant__ uint16_T gC_pos_X[MAX_IND_READ];
__constant__ uint16_T gC_pos_Y[MAX_IND_READ];


#define MAX(x,y) (x>y?x:y);
#define MIN(x,y) (x<y?x:y);
#define ABS(x) (x>0?x:-x);


int checkLastError(char * msg)
{
    cudaError_t cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        char err[512];
        sprintf(err, "getprojection failed \n %s: %s. \n", msg, cudaGetErrorString(cudaStatus));
        mexPrintf(err);
        return 1;
    }
    return 0;
}


/*
 * Device code
 */
/*********** fast const memory based version ***************/ 

__global__ void readFromArray_c_fast(float2 * sarray, const float2 * larray, 
                   cuint Np_px,cuint Np_py, cuint Np_pz,cuint Np_ox, cuint Np_oy, 
                                        cuint Npos)  {
    // Location in a 3D matrix
    int idx= blockIdx.x * blockDim.x + threadIdx.x;
    int idy= blockIdx.y * blockDim.y + threadIdx.y;
    int id = blockIdx.z * blockDim.z + threadIdx.z;

    if ( idx < Np_px & idy < Np_py & id < Npos)
    {
        int idz = gC_ind_read[id]-1;  // go only through some of the indices 
        int id_large =    gC_pos_X[idz]+idx + Np_ox*(gC_pos_Y[idz]+idy);
        int id_small = idx + Np_px*idy + Np_px*Np_py*idz ;
        sarray[ id_small ].x =  larray[ id_large ].x ;
        sarray[ id_small ].y = larray[ id_large ].y ;
    }
}
/*********** global memory based version ***************/ 

__global__ void readFromArray_c(float2 * sarray, const float2 * larray, cuint16* ind_read, cuint16* pos_X, cuint16* posY, 
                   cuint Np_px,cuint Np_py, cuint Np_pz,cuint Np_ox, cuint Np_oy, 
                                        cuint Npos)  {
    // Location in a 3D matrix
    int idx= blockIdx.x * blockDim.x + threadIdx.x;
    int idy= blockIdx.y * blockDim.y + threadIdx.y;
    int id = blockIdx.z * blockDim.z + threadIdx.z;

    if ( idx < Np_px & idy < Np_py & id < Npos)
    {
        int idz = ind_read[id]-1;  // go only through some of the indices 
        int id_large =    pos_X[idz]+idx + Np_ox*(posY[idz]+idy);
        int id_small = idx + Np_px*idy + Np_px*Np_py*idz ;
        sarray[ id_small ].x =  larray[ id_large ].x ;
        sarray[ id_small ].y = larray[ id_large ].y ;
    }
}


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    int i ; 
    char const * const errId = "parallel:gpu:mexGPUExample:InvalidInput";
    char const * const errMsg = "Invalid input to MEX file.";

    /* Initialize the MathWorks GPU API. */
    //mxInitGPU();
    
    /* Check for proper number of arguments. */
    if (nrhs != 5) 
        mexErrMsgTxt("Five input arguments required");
  

  for (i=0; i < 2; i++) {
      if ( !mxIsGPUArray(prhs[i]) && !mxIsCell(prhs[i]) ){
          printf("Input %d is not GPU array / cell \n",i+1);
           mexErrMsgIdAndTxt("MexError:ptycho","Inputs must be of correct type.");
      }
  }


 // load positions 

 
     const mxGPUArray * m_positions_x = mxGPUCreateFromMxArray(prhs[2]); 
     if ((mxGPUGetClassID(m_positions_x) != mxUINT16_CLASS)) {
           mexPrintf("m_positions_x\n");
           mexErrMsgIdAndTxt(errId, errMsg);
     }
     const cuint16 * p_positions_x = (cuint16 *)mxGPUGetDataReadOnly(m_positions_x);

     const mxGPUArray * m_positions_y = mxGPUCreateFromMxArray(prhs[3]); 
     if ((mxGPUGetClassID(m_positions_y) != mxUINT16_CLASS)) {
           mexPrintf("m_positions_y\n");
           mexErrMsgIdAndTxt(errId, errMsg);
     }
     const cuint16 * p_positions_y = (cuint16 *)mxGPUGetDataReadOnly(m_positions_y);

  
    /**** copy of the array is the slowest operation *****/ 
    // Now it is writting directly into the input field !!! 
    //mxGPUArray * m_obj_proj = mxGPUCopyFromMxArray(prhs[0]); 
    mxGPUArray * m_obj_proj = const_cast<mxGPUArray*>(mxGPUCreateFromMxArray(prhs[0]));  // 
    if ((mxGPUGetClassID(m_obj_proj) != mxSINGLE_CLASS)) {
           mexPrintf("m_obj_proj\n");
           mexErrMsgIdAndTxt(errId, errMsg);
     }
     if (mxGPUGetComplexity(m_obj_proj) != mxCOMPLEX) {
        mexPrintf("m_obj_proj is not complex \n");
           mexErrMsgIdAndTxt(errId, errMsg);
     }
     float2 * p_obj_proj = (float2 *)mxGPUGetData(m_obj_proj);


    if (!mxIsCell(prhs[1]) || !mxIsCell(prhs[4]))
           mexErrMsgIdAndTxt("MexError:ptycho","Object and indices has to be in cell/cells !! ");
    if (mxGetNumberOfElements(prhs[1]) != mxGetNumberOfElements(prhs[4]))
           mexErrMsgIdAndTxt("MexError:ptycho","Number of objects != number of indices ");


    /* Get dimension of probe and object */
    const unsigned int  Ndims = (unsigned int)mxGPUGetNumberOfDimensions(m_obj_proj);


    cuint Ncells = mxGetNumberOfElements(prhs[1]);
    const unsigned int Np_pp = mxGPUGetNumberOfElements(m_positions_x);
    

    if (Np_pp < MAX_IND_READ) {
            cudaMemcpyToSymbol(gC_pos_X, p_positions_x, Np_pp*sizeof(uint16_T), 0, cudaMemcpyHostToDevice);
            cudaMemcpyToSymbol(gC_pos_Y, p_positions_y, Np_pp*sizeof(uint16_T), 0, cudaMemcpyHostToDevice);
            checkLastError("after cudaMemcpyToSymbol pos");
        } 

    for (int l=0; l<Ncells; l++)
    {
        // read the cell content 
        mxArray * mx_object = mxGetCell(prhs[1],l);
        mxArray * mx_ind = mxGetCell(prhs[4],l);
        int N_ok = mxGetNumberOfElements(mx_ind);
        if(N_ok == 0)
            continue;

        const  mxGPUArray * m_object = mxGPUCreateFromMxArray(mx_object); 
         if ((mxGPUGetClassID(m_object) != mxSINGLE_CLASS)) {
               mexPrintf("m_object\n");
               mexErrMsgIdAndTxt(errId, errMsg);
         }
        if (mxGPUGetComplexity(m_object) != mxCOMPLEX) {
            mexPrintf("m_object is not complex \n");
               mexErrMsgIdAndTxt(errId, errMsg);
         }
         const float2 * p_object = (float2 *)mxGPUGetDataReadOnly(m_object);


         const mxGPUArray * m_ind_ok = mxGPUCreateFromMxArray(mx_ind); 
         if ((mxGPUGetClassID(m_ind_ok) != mxUINT16_CLASS)) {
               mexPrintf("m_ind_ok\n");
               mexErrMsgIdAndTxt(errId, errMsg);
         }
         const cuint16 *  p_ind_ok = (cuint16 *)mxGPUGetDataReadOnly(m_ind_ok);

         const mwSize * Np_o = mxGPUGetDimensions(m_object);
         const mwSize * Np_p = mxGPUGetDimensions(m_obj_proj);
         const unsigned int Npos = mxGPUGetNumberOfElements(m_ind_ok);



//       mexPrintf("Ndims %i  Np_o %i %i  Np_p %i %i %i  Npos %i \n " ,Ndims,Np_o[0],Np_o[1],Np_p[0],Np_p[1],Np_p[2],Npos);




        if (Ndims == 3 && Npos > Np_p[2]) {
            printf("wrong size of update / positions  %i", Ndims);
            mexErrMsgIdAndTxt("MexError:ptycho","wrong size of update / positions");
        } 



        // Choose a reasonably sized number of threads in each dimension for the block.
        int const threadsPerBlockEachDim =  32;
        // Compute the thread block and grid sizes based on the board dimensions.
        int const blocksPerGrid_M = (Np_p[0] + threadsPerBlockEachDim - 1) / threadsPerBlockEachDim;
        int const blocksPerGrid_N = (Np_p[1] + threadsPerBlockEachDim - 1) / threadsPerBlockEachDim;
        int const blocksPerGrid_O = Npos;

        //  mexPrintf("Threads %i %i %i \n ", blocksPerGrid_M, blocksPerGrid_N, blocksPerGrid_O);

        dim3 const dimBlock(blocksPerGrid_M, blocksPerGrid_N, blocksPerGrid_O);
        dim3 const dimThread(threadsPerBlockEachDim, threadsPerBlockEachDim, 1);


        checkLastError("after dimThread");

        //mexPrintf("Blocks %i %i %i \n ", dimThread.x, dimThread.y, dimThread.z);

        if (Np_pp > MAX_IND_READ) {
            //mexPrintf( "More than %i positions may be slow \n", MAX_IND_READ);
        } else {
            cudaMemcpyToSymbol(gC_ind_read, p_ind_ok, Npos*sizeof(uint16_T), 0, cudaMemcpyHostToDevice);
            checkLastError("after cudaMemcpyToSymbol pos");
        } 

        checkLastError("after cudaMemcpyToSymbol");


        //=============  run the kernel ======================
        if (Npos < MAX_IND_READ) 
            readFromArray_c_fast<<<dimBlock, dimThread>>>(p_obj_proj,p_object, Np_p[0],Np_p[1],Np_p[2],Np_o[0],Np_o[1], Npos);
        else
            readFromArray_c<<<dimBlock, dimThread>>>(p_obj_proj,p_object,p_ind_ok, p_positions_x,p_positions_y,Np_p[0],Np_p[1],Np_p[2],Np_o[0],Np_o[1], Npos);

    
        checkLastError("after kernel");


        mxGPUDestroyGPUArray(m_object);
        mxGPUDestroyGPUArray(m_ind_ok);


    }

    cudaThreadSynchronize();


    // plhs[0] = mxGPUCreateMxArrayOnGPU(m_obj_proj);
    mxGPUDestroyGPUArray(m_obj_proj);
    mxGPUDestroyGPUArray(m_positions_x);
    mxGPUDestroyGPUArray(m_positions_y);

    
  return;
}




