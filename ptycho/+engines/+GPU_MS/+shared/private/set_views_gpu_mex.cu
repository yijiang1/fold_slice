/*  
    Set complex views to complex object 
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
__constant__ uint16_T gC_ind_read[MAX_IND_READ_DEV];
__constant__ uint16_T gC_pos_X[MAX_IND_READ_DEV];
__constant__ uint16_T gC_pos_Y[MAX_IND_READ_DEV];


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

/***********  reduction of object projection array to single object ***************/ 
template <bool useGlobal>
__global__ void addToArray_r( float const * sarray, float * larray,  cuint16* ind_read, cuint16* pos_X, cuint16* posY, 
                            cuint Np_px,cuint Np_py, cuint Np_pz,cuint Np_ox, cuint Np_oy, 
                                        cuint Npos, const bool isFlat)  {
    // Location in a 3D matrix
    int idx= blockIdx.x * blockDim.x + threadIdx.x;
    int idy= blockIdx.y * blockDim.y + threadIdx.y;
    if ( idx < Np_px & idy < Np_py ) {
        int idz, id_large; 
        float sarray_val; 
        for(int id = 0; id < Npos; id++)
        {   
            if (useGlobal) {
                //  fast const memory based version
                idz = gC_ind_read[id]-1;  // go only through some of the indices 
                id_large =  gC_pos_X[idz]+idx + Np_ox*(gC_pos_Y[idz]+idy);
            } else {
                // slower global memory based version
                idz = ind_read[id]-1;  // go only through some of the indices 
                id_large =  pos_X[idz]+idx + Np_ox*(posY[idz]+idy);
            }
            int id_small = idx + Np_px*idy ;
            if (!isFlat)
                id_small = id_small + Np_px*Np_py*idz ;
            
            // prevent extra memory load 
            sarray_val = (isFlat && (idz > 0)) ? sarray_val: sarray[ id_small ];

            //larray[id_large]  += sarray_val;
            //__syncthreads(); 

            // slowest step, without atomicAdd it misses some values  
            atomicAdd(&larray[id_large] ,sarray_val);
        }
    }
}

template <bool useGlobal>
__global__ void addToArray_c( float2 const * sarray, float2 * larray,  cuint16* ind_read, cuint16* pos_X, cuint16* posY, 
                            cuint Np_px,cuint Np_py, cuint Np_pz,cuint Np_ox, cuint Np_oy, 
                                        cuint Npos, const bool isFlat)  {
    // Location in a 3D matrix
    int idx= blockIdx.x * blockDim.x + threadIdx.x;
    int idy= blockIdx.y * blockDim.y + threadIdx.y;
    if ( idx < Np_px & idy < Np_py ) {
        int idz, id_large; 
        float2 sarray_val; 
        for(int id = 0; id < Npos; id++)
        {   
            if (useGlobal) {
                //  fast const memory based version
                idz = gC_ind_read[id]-1;  // go only through some of the indices 
                id_large =  gC_pos_X[idz]+idx + Np_ox*(gC_pos_Y[idz]+idy);
            } else {
                // slower global memory based version
                idz = ind_read[id]-1;  // go only through some of the indices 
                id_large =  pos_X[idz]+idx + Np_ox*(posY[idz]+idy);
            }
            int id_small = idx + Np_px*idy ;
            if (!isFlat)
                id_small = id_small + Np_px*Np_py*idz ;

            // prevent extra memory load 
            sarray_val = (isFlat && (idz > 0)) ? sarray_val: sarray[ id_small ];


            //larray[id_large].x  += sarray_val.x;
            //larray[id_large].y  += sarray_val.y;
            //__syncthreads(); 

            // slowest step, without atomicAdd it misses some values  
            atomicAdd(&larray[id_large].x ,sarray_val.x);
            atomicAdd(&larray[id_large].y ,sarray_val.y);
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
    if (nrhs != 5) 
    mexErrMsgTxt("Five input arguments required");
  

    // Input must be of type single. 
  for (i=0; i < 2; i++) {
      if ( !mxIsGPUArray(prhs[i]) && !mxIsCell(prhs[i]) ){
          printf("Input %d is not cell array\n",i+1);
           mexErrMsgIdAndTxt("MexError:ptycho","Inputs must be of correct type.");
      }
  }
    // Input must be of type int16.
  for (i=2; i<4; i++){
    if (mxIsUint16(prhs[i]) != 1){
          printf("Input %d is not integer\n",i+1);
          mexErrMsgIdAndTxt("MexError:ptycho","Inputs must be of correct type uint16.");
    }
  }
 
   // It cannot be one-dimensional 
  if(mxGetNumberOfDimensions(prhs[0]) < 2) {
    printf("The 1st input argument must have at least two dimensions.");
    mexErrMsgIdAndTxt("MexError:ptycho","wrong number of dimensions");
  }
    // It cannot be more than 3-dimensional 
   if(mxGetNumberOfDimensions(prhs[0]) > 3) {
    printf("The 1st input argument must have at most three dimensions.");
    mexErrMsgIdAndTxt("MexError:ptycho","wrong number of dimensions");
  }
  // Check that arrays are complex

     const mxGPUArray * m_obj_proj = mxGPUCreateFromMxArray(prhs[0]); 
     if ((mxGPUGetClassID(m_obj_proj) != mxSINGLE_CLASS)) {
           mexPrintf("m_obj_proj\n");
           mexErrMsgIdAndTxt(errId, errMsg);
     }

     const mxGPUArray * m_positions_x = mxGPUCreateFromMxArray(prhs[2]); 
     if ((mxGPUGetClassID(m_positions_x) != mxUINT16_CLASS)) {
           mexPrintf("m_positions_x\n");
           mexErrMsgIdAndTxt(errId, errMsg);
     }
     const uint16_T * p_positions_x = (uint16_T *)mxGPUGetDataReadOnly(m_positions_x);

     const mxGPUArray * m_positions_y = mxGPUCreateFromMxArray(prhs[3]); 
     if ((mxGPUGetClassID(m_positions_y) != mxUINT16_CLASS)) {
           mexPrintf("m_positions_y\n");
           mexErrMsgIdAndTxt(errId, errMsg);
     }
     const uint16_T * p_positions_y = (uint16_T *)mxGPUGetDataReadOnly(m_positions_y);

    if (!mxIsCell(prhs[1]) || !mxIsCell(prhs[4]))
           mexErrMsgIdAndTxt("MexError:ptycho","Object and indices has to be in cell/cells !! ");
    if (mxGetNumberOfElements(prhs[1]) != mxGetNumberOfElements(prhs[4]))
           mexErrMsgIdAndTxt("MexError:ptycho","Number of objects != number of indices ");

    cuint Ncells = mxGetNumberOfElements(prhs[1]);
    const unsigned int Np_pp = mxGPUGetNumberOfElements(m_positions_y);

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

         mxGPUArray * m_object = const_cast<mxGPUArray*>(mxGPUCreateFromMxArray(mx_object)); //    mxGPUCopyFromMxArray(prhs[1]); 
         if ((mxGPUGetClassID(m_object) != mxSINGLE_CLASS)) {
               mexPrintf("m_object\n");
               mexErrMsgIdAndTxt(errId, errMsg);
         }

        if (mxGPUGetComplexity(m_object) != mxGPUGetComplexity(m_obj_proj)) {
            mexPrintf("m_object and m_obj_proj complexity has to be the same \n");
               mexErrMsgIdAndTxt(errId, errMsg);
         }



         const mxGPUArray * m_ind_ok = mxGPUCreateFromMxArray(mx_ind); 
         if ((mxGPUGetClassID(m_ind_ok) != mxUINT16_CLASS)) {
               mexPrintf("m_ind_ok\n");
               mexErrMsgIdAndTxt(errId, errMsg);
         }
         const uint16_T *  p_ind_ok = (uint16_T *)mxGPUGetDataReadOnly(m_ind_ok);



        // Get dimension of probe and object 
        const unsigned int  Ndims = (unsigned int)mxGPUGetNumberOfDimensions(m_obj_proj);
        const mwSize * Np_o = mxGPUGetDimensions(m_object);
        const mwSize * Np_p = mxGPUGetDimensions(m_obj_proj);
        const unsigned int Npos = mxGPUGetNumberOfElements(m_ind_ok);

        //mexPrintf("Ndims %i  Np_o %i %i  Np_p %i %i %i  Npos %i \n " ,Ndims,Np_o[0],Np_o[1],Np_p[0],Np_p[1],Np_p[3],Npos);




      if (Ndims == 3 && Npos > Np_p[2]) {
          printf("wrong size of update / positions  %i", Ndims);
          mexErrMsgIdAndTxt("MexError:ptycho","wrong size of update / positions");
      } 


        if (Npos > MAX_IND_READ) {
            //mexPrintf( "More than %i positions may be slow \n", MAX_IND_READ);
        } else {
            cudaMemcpyToSymbol(gC_ind_read, p_ind_ok, Npos*sizeof(uint16_T), 0, cudaMemcpyHostToDevice);
            checkLastError("after cudaMemcpyToSymbol pos");
        }


        // Choose a reasonably sized number of threads in each dimension for the block.
        int const threadsPerBlockEachDim =  32;
        // Compute the thread block and grid sizes based on the board dimensions.
        int const blocksPerGrid_M = (Np_p[0] + threadsPerBlockEachDim - 1) / threadsPerBlockEachDim;
        int const blocksPerGrid_N = (Np_p[1] + threadsPerBlockEachDim - 1) / threadsPerBlockEachDim;
        int const blocksPerGrid_O = 1;

        dim3 const dimBlock(blocksPerGrid_M, blocksPerGrid_N, blocksPerGrid_O);
        dim3 const dimThread(threadsPerBlockEachDim, threadsPerBlockEachDim, 1);

        checkLastError("after dimThread");

        checkLastError("after cudaMemcpyToSymbol");

        const bool isFlat = (Ndims == 2);
        const bool isComplex = mxGPUGetComplexity(m_obj_proj) == mxCOMPLEX;
        // ================== call the right kernel =================== 
        if (isComplex) {
            const float2 * p_obj_proj = (float2 *)mxGPUGetDataReadOnly(m_obj_proj);
            float2 * p_object = (float2 *)mxGPUGetData(m_object);
            if (Np_pp < MAX_IND_READ) 
                addToArray_c<true><<<dimBlock, dimThread>>>(p_obj_proj,p_object, p_ind_ok, p_positions_x,p_positions_y ,Np_p[0],Np_p[1],Np_p[2],Np_o[0],Np_o[1], Npos, isFlat);
            else
                addToArray_c<false><<<dimBlock, dimThread>>>(p_obj_proj,p_object, p_ind_ok, p_positions_x,p_positions_y ,Np_p[0],Np_p[1],Np_p[2],Np_o[0],Np_o[1], Npos, isFlat);
        } else {
             const float * p_obj_proj = (float *)mxGPUGetDataReadOnly(m_obj_proj);
             float * p_object = (float *)mxGPUGetData(m_object);
            if (Np_pp < MAX_IND_READ) 
                addToArray_r<true><<<dimBlock, dimThread>>>(p_obj_proj,p_object, p_ind_ok, p_positions_x,p_positions_y ,Np_p[0],Np_p[1],Np_p[2],Np_o[0],Np_o[1], Npos, isFlat);
            else
                addToArray_r<false><<<dimBlock, dimThread>>>(p_obj_proj,p_object, p_ind_ok, p_positions_x,p_positions_y ,Np_p[0],Np_p[1],Np_p[2],Np_o[0],Np_o[1], Npos, isFlat);
        }

        checkLastError("after kernel");

        mxGPUDestroyGPUArray(m_object);
        mxGPUDestroyGPUArray(m_ind_ok);
    }

    cudaThreadSynchronize();

    // plhs[0] = mxGPUCreateMxArrayOnGPU(m_object);
    mxGPUDestroyGPUArray(m_obj_proj);
    mxGPUDestroyGPUArray(m_positions_x);
    mxGPUDestroyGPUArray(m_positions_y);

    
  return;
}




