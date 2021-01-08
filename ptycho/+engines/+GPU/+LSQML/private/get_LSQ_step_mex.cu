/*  
    Set complex views to complex object 
    mexcuda -output +engines/+GPU/get_optimal_LSQ_step_mex +engines/+GPU/get_optimal_LSQ_step_mex.cu
 */

#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <math.h> 
#include <stdio.h> 
#include <iostream>
#include <list>

typedef const unsigned int   cuint;
typedef const uint16_T cuint16;
#define MAX_BLOCK_DIM_SIZE 65535


/*
 * Device code
 */

// allocate shared memory so that all functions can see it 
extern __shared__ float sdata[];
const unsigned int MAX_IND_READ = 10000;
__constant__ uint8_T gC_pind[MAX_IND_READ];


int checkLastError(char * msg)
{
    cudaError_t cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        char err[512];
        sprintf(err, "get_optimal_LSQ_step_ker failed \n %s: %s. \n", msg, cudaGetErrorString(cudaStatus));
        mexPrintf(err);
        return 1;
    }
    return 0;
}


/*********** fast inplace version of LSQ step calculation *************/

template <unsigned int blockSize>
__device__ void calculate_AA_matrix( const float2 *P_f, const float2 *O_f,const float2 *dP_f, const float2 *dO_f, const float2 *chi_f, 
                                   float &AA1, float2 &AA2,float2 &AA3, float &AA4, 
                                   float  &Atb1, float &Atb2 , const float lambda, 
                                    cuint Np_x, cuint Np_y, cuint Npixz, cuint idz, 
                                    cuint Nblocks, const bool single_probe, cuint id, cuint tid)
{
    
    float2 dO, dP, O, P, chi; 
    // load to local memory 
    cuint id3 = id + Np_x*Np_y*idz;
    O   = O_f[id3]  ;
    dO  = dO_f[id3] ; 
    chi = chi_f[id3];

    if (single_probe) {
        // single shared 2D probe 
        P  = P_f[id];
        dP = dP_f[id];
     } else {
        // unshared probe => size(dp,3) == Nscans 
        dP = dP_f[id + (gC_pind[idz]-1)*Np_x*Np_y]; 
        // position in 3D array 
        P  = P_f[id3];   
    }
    
    // make auxiliary variables 
    float2 dOP, dPO, cdPO, cdOP; 

    //  dOP = dO.*P;
    dOP.x = dO.x * P.x -  dO.y * P.y; 
    dOP.y = dO.y * P.x +  dO.x * P.y;

    //  dPO = dP.*O;
    dPO.x = dP.x * O.x -  dP.y * O.y; 
    dPO.y = dP.y * O.x +  dP.x * O.y;

    // cdOP = conj(dOP);
    cdOP.x =  dOP.x; 
    cdOP.y = -dOP.y; 

    // cdPO = conj(dPO);
    cdPO.x =  dPO.x; 
    cdPO.y = -dPO.y; 

    // AA1 = abs(dOP).^2+lambda;
    AA1 = dOP.x * dOP.x + dOP.y * dOP.y ;

    // AA2 = (dOP .* cdPO);
    AA2.x = dOP.x * cdPO.x - dOP.y * cdPO.y ;
    AA2.y = dOP.x * cdPO.y + dOP.x * cdPO.y ; 

    // AA3 = conj(AA2);
    AA3.x = AA2.x; 
    AA3.y = -AA2.y; 

    // AA4 = abs(dPO)^2+lambda;
    AA4 = dPO.x * dPO.x + dPO.y * dPO.y ;

    // Atb1 = real(cdOP .* chi);
    Atb1 =  cdOP.x*chi.x - cdOP.y*chi.y; 

    // Atb2 = real(cdPO .* chi);
    Atb2 = cdPO.x*chi.x -  cdPO.y*chi.y; 


    // add to the shared gpu memory 
    sdata[tid            ] = AA1;
    sdata[tid+  blockSize] = AA2.x;
    sdata[tid+2*blockSize] = AA2.y;
    sdata[tid+3*blockSize] = AA3.x;
    sdata[tid+4*blockSize] = AA3.y;
    sdata[tid+5*blockSize] = AA4;
    sdata[tid+6*blockSize] = Atb1;
    sdata[tid+7*blockSize] = Atb2;
}

template <unsigned int blockSize>
__device__ void add_to_shared_array( cuint tid, cuint offset )
{
    // another loop unrolling 
    sdata[tid + 0*blockSize] += sdata[tid + 0*blockSize + offset];
    sdata[tid + 1*blockSize] += sdata[tid + 1*blockSize + offset];
    sdata[tid + 2*blockSize] += sdata[tid + 2*blockSize + offset];
    sdata[tid + 3*blockSize] += sdata[tid + 3*blockSize + offset];
    sdata[tid + 4*blockSize] += sdata[tid + 4*blockSize + offset];
    sdata[tid + 5*blockSize] += sdata[tid + 5*blockSize + offset];
    sdata[tid + 6*blockSize] += sdata[tid + 6*blockSize + offset];
    sdata[tid + 7*blockSize] += sdata[tid + 7*blockSize + offset];
}


template <unsigned int blockSize>
__device__ void reduce_shared_array( cuint tid )
{
    // do reduction in shared mem using unrolled loops 
    if (blockSize >= 1024){ if (tid < 512) { add_to_shared_array<blockSize>(tid,512); } __syncthreads(); }
    if (blockSize >= 512) { if (tid < 256) { add_to_shared_array<blockSize>(tid,256); } __syncthreads(); }
    if (blockSize >= 256) { if (tid < 128) { add_to_shared_array<blockSize>(tid,128); } __syncthreads(); }
    if (blockSize >= 128) { if (tid <  64) { add_to_shared_array<blockSize>(tid,64); } __syncthreads(); }
    // why not do the same for all 
    if (blockSize >= 64) { if (tid <  32) { add_to_shared_array<blockSize>(tid,32); } __syncthreads(); }
    if (blockSize >= 32) { if (tid <  16) { add_to_shared_array<blockSize>(tid,16); } __syncthreads(); }
    if (blockSize >= 16) { if (tid <  8)  { add_to_shared_array<blockSize>(tid,8); } __syncthreads(); }
    if (blockSize >= 8)  { if (tid <  4)  { add_to_shared_array<blockSize>(tid,4); } __syncthreads(); }
    if (blockSize >= 4)  { if (tid <  2)  { add_to_shared_array<blockSize>(tid,2); } __syncthreads(); }
    if (blockSize >= 2)  { if (tid <  1)  { add_to_shared_array<blockSize>(tid,1); } __syncthreads(); }
}


// fast kernel for estimation of optimal probe and object steps 

template <unsigned int blockSize>
__global__ void get_optimal_LSQ_step_ker( float2 const * P_f, float2 const * O_f,float2 const * dP_f, float2 const * dO_f,
                            float2 const * chi_f, const float lambda, 
                            float2 * AA, float * Atb, cuint Np_x,cuint Np_y, cuint Npixz, cuint Nblocks, const bool single_probe)  {

const mwSize tid = threadIdx.x;
// do only every second block 
//cuint i = blockIdx.x*(blockSize*2) + threadIdx.x;
const mwSize i = blockIdx.x*(blockDim.x) + threadIdx.x;

const mwSize N2 = Np_x*Np_y;
mwSize AA_page_id, Atb_page_id;

float2 AA2, AA3;
float AA1, AA4, Atb1, Atb2;

for(int n = 0; n < 8; n++)
    sdata[tid + n*blockSize ] = 0 ;


if(i < N2)
{

    //for(int n = 0; n < 8*blockSize; n++)
    //    sdata[n] = 0 ;


    // Page in a 3D matrix
    for(int idz = 0; idz < Npixz; idz++)
    {
        unsigned int ii = i ; 

        // empty the share memory 
        for(int n = 0; n < 8; n++)
            sdata[tid + n*blockSize ] = 0 ;



        // get coeficients for the AA matrix + right size Atb vector and add them to the shared array 
        calculate_AA_matrix<blockSize>(P_f ,O_f ,dP_f ,dO_f ,chi_f , 
                    AA1, AA2, AA3, AA4, Atb1, Atb2, lambda, Np_x, Np_y, Npixz, idz, Nblocks,single_probe, ii, tid); 


        __syncthreads();

   
        // reduce the shared memory data 
        reduce_shared_array<blockSize>( tid ); 
   

        // write result for this block to global mem 
        if (tid == 0)     {
            // store data to the AA matrix 
            AA_page_id = 4*idz + 4*Npixz*blockIdx.x; 
            Atb_page_id = 2*idz + 2*Npixz*blockIdx.x; 

            AA[ 0 + AA_page_id].x = sdata[0*blockSize];
            AA[ 0 + AA_page_id].y = 0;   // needs to be set to zero or initalized to zero when created 

            AA[ 1 + AA_page_id].x = sdata[1*blockSize];
            AA[ 1 + AA_page_id].y = sdata[2*blockSize]; 
            AA[ 2 + AA_page_id].x = sdata[3*blockSize]; 
            AA[ 2 + AA_page_id].y = sdata[4*blockSize];
            AA[ 3 + AA_page_id].x = sdata[5*blockSize];
            AA[ 3 + AA_page_id].y = 0;  

            Atb[ 0 + Atb_page_id] = sdata[6*blockSize]; 
            Atb[ 1 + Atb_page_id] = sdata[7*blockSize];  

        }
    }
}
}





unsigned int nextPow2( unsigned int x ) {
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
}

void getNumBlocksAndThreads(int n, int maxBlocks, int maxThreads, int &blocks, int &threads)
{
    threads = (n < maxThreads*2) ? nextPow2((n + 1)/ 2) : maxThreads;
    blocks = (n + (threads * 2 - 1)) / (threads * 2);
    blocks = min(maxBlocks, blocks);
}



void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    char const * const errId = "parallel:gpu:mexGPUExample:InvalidInput";
    char const * const errMsg = "Invalid input to MEX file.";

    
    // Check for proper number of arguments. 

    if (nrhs != 7) 
    mexErrMsgTxt("Seven input arguments required");



     const mxGPUArray * m_chi = mxGPUCreateFromMxArray(prhs[0]); 
     if ((mxGPUGetClassID(m_chi) != mxSINGLE_CLASS) || (mxGPUGetComplexity(m_chi) != mxCOMPLEX)) {
           mexPrintf("m_chi\n");
           mexErrMsgIdAndTxt(errId, errMsg);
     }
     const float2 * p_chi = (float2 *)mxGPUGetDataReadOnly(m_chi);

     const mxGPUArray * m_dO = mxGPUCreateFromMxArray(prhs[1]); 
     if ((mxGPUGetClassID(m_dO) != mxSINGLE_CLASS) || (mxGPUGetComplexity(m_dO) != mxCOMPLEX)) {
           mexPrintf("m_dO\n");
           mexErrMsgIdAndTxt(errId, errMsg);
     }
     const float2 * p_dO = (float2 *)mxGPUGetDataReadOnly(m_dO);

     const mxGPUArray * m_dP = mxGPUCreateFromMxArray(prhs[2]); 
     if ((mxGPUGetClassID(m_dP) != mxSINGLE_CLASS) || (mxGPUGetComplexity(m_dP) != mxCOMPLEX)) {
           mexPrintf("m_dP\n");
           mexErrMsgIdAndTxt(errId, errMsg);
     }
     const float2 * p_dP = (float2 *)mxGPUGetDataReadOnly(m_dP);

     const mxGPUArray * m_O = mxGPUCreateFromMxArray(prhs[3]); 
     if ((mxGPUGetClassID(m_O) != mxSINGLE_CLASS) || (mxGPUGetComplexity(m_O) != mxCOMPLEX)) {
           mexPrintf("m_O\n");
           mexErrMsgIdAndTxt(errId, errMsg);
     }
     const float2 * p_O = (float2 *)mxGPUGetDataReadOnly(m_O);

     const mxGPUArray * m_P = mxGPUCreateFromMxArray(prhs[4]); 
     if ((mxGPUGetClassID(m_P) != mxSINGLE_CLASS) || (mxGPUGetComplexity(m_P) != mxCOMPLEX)) {
           mexPrintf("m_P\n");
           mexErrMsgIdAndTxt(errId, errMsg);
     }
     const float2 * p_P = (float2 *)mxGPUGetDataReadOnly(m_P);

     const mxGPUArray * m_P_ind = mxGPUCreateFromMxArray(prhs[6]); 
     if (mxGPUGetClassID(m_P_ind) != mxUINT8_CLASS)  {
           mexPrintf("m_P_ind class %i\n", mxGPUGetClassID(m_P_ind));
           mexErrMsgIdAndTxt(errId, errMsg);
     }
     const uint8_T * p_P_ind = (uint8_T *)mxGPUGetDataReadOnly(m_P_ind);
     const unsigned int Npos = mxGPUGetNumberOfElements(m_P_ind);

    if (Npos > MAX_IND_READ) {
       mexErrMsgIdAndTxt(errId, "Maximal size of input block exceeded");
    }


    // Get dimension of probe and object 
    const unsigned int  Ndims = (unsigned int)mxGPUGetNumberOfDimensions(m_chi);
    if (Ndims != 3) {
        mexErrMsgIdAndTxt(errId, "Inputs has to be 3 dimensional\n");
     }
    const mwSize * Npix = mxGPUGetDimensions(m_chi);
    const mwSize * Npix_probe = mxGPUGetDimensions(m_P);
    const mwSize * Npix_probe_upd = mxGPUGetDimensions(m_dP);
    const mwSize  Ndims_probe = mxGPUGetNumberOfDimensions(m_P);
    const mwSize  Ndims_probe_upd = mxGPUGetNumberOfDimensions(m_dP);

    if ((Npix[2] != Npos)) {
       mexErrMsgIdAndTxt(errId, "Number of probe indices has to match size of inputs (%i vs %i) \n", Npix[2], Npos);
    }
    if ((Npix_probe[2] != Npix[2]) && (Ndims_probe != 2)) {
       mexErrMsgIdAndTxt(errId, "Dimension of probe has to match size of inputs (%i vs %i) \n", Npix_probe[2], Npix[2]);
    }


    float lambda = mxGetScalar(prhs[5]);

    
    
    const bool single_probe =Ndims_probe == 2  ;

   

    cudaMemcpyToSymbol(gC_pind, p_P_ind, Npos*sizeof(uint8_T), 0, cudaMemcpyHostToDevice);
    checkLastError("after cudaMemcpyToSymbol pos");

    // Choose a reasonably sized number of threads in each dimension for the block.
    int maxThreads = 1024;  // number of threads per block, does not work with 1024, I dont know why
    int threads = 0, blocks = 0;

    cuint n = Npix[0]*Npix[1]; 

    threads = (n < maxThreads) ? nextPow2((n + 1)/ 2) : maxThreads;
    blocks = (n + (threads  - 1)) / (threads );



    dim3 dimBlock(threads, 1, 1);
    dim3 dimGrid(blocks, 1, 1);
    
    // allocation size needed for shared GPU memory , it needs to reduce 2 float2 elements and 4 float elements
    int smemSize = 8* threads * sizeof(float);
    
    //mexPrintf("threads %i blocks %i smemSize %i \n", threads, blocks, smemSize);


	// allocate output fields  
    mwSize matrix_size[4] = {2,2,Npix[2],blocks};
    mxGPUArray * m_AA = mxGPUCreateGPUArray(
        4,
        matrix_size,
        mxSINGLE_CLASS,
        mxCOMPLEX,
        MX_GPU_DO_NOT_INITIALIZE);  // MX_GPU_DO_NOT_INITIALIZE , MX_GPU_INITIALIZE_VALUES
    float2 * p_AA = (float2 *)mxGPUGetData(m_AA);

    mwSize vector_size[4] = {2,1,Npix[2],blocks};
    mxGPUArray * m_Atb = mxGPUCreateGPUArray(
        4,
        vector_size,
        mxSINGLE_CLASS,
        mxREAL,
        MX_GPU_DO_NOT_INITIALIZE);
    float * p_Atb = (float *)mxGPUGetData(m_Atb);


  
    checkLastError("after dimThread");
    switch (threads)
    {
        case 1024:
           get_optimal_LSQ_step_ker< 1024><<< dimGrid, dimBlock, smemSize>>>(  p_P, p_O,p_dP, p_dO, p_chi, lambda, 
                                                    p_AA, p_Atb, Npix[0], Npix[1], Npix[2], blocks, single_probe);
           break;
        case 512:
           get_optimal_LSQ_step_ker< 512><<< dimGrid, dimBlock, smemSize>>>(  p_P, p_O,p_dP, p_dO, p_chi, lambda, 
                                                    p_AA, p_Atb, Npix[0], Npix[1], Npix[2], blocks, single_probe);
           break;
        case 256:
           get_optimal_LSQ_step_ker< 256><<< dimGrid, dimBlock, smemSize>>>(  p_P, p_O,p_dP, p_dO, p_chi, lambda, 
                                                    p_AA, p_Atb, Npix[0], Npix[1], Npix[2], blocks, single_probe);
            break;
        case 128:
           get_optimal_LSQ_step_ker< 128><<< dimGrid, dimBlock, smemSize>>>(  p_P, p_O,p_dP, p_dO, p_chi, lambda, 
                                                    p_AA, p_Atb, Npix[0], Npix[1], Npix[2], blocks, single_probe);
            break;
        case 64:
           get_optimal_LSQ_step_ker< 64><<< dimGrid, dimBlock, smemSize>>>(  p_P, p_O,p_dP, p_dO, p_chi, lambda, 
                                                   p_AA, p_Atb, Npix[0], Npix[1], Npix[2], blocks, single_probe);
            break;
        case 32:
           get_optimal_LSQ_step_ker< 32><<< dimGrid, dimBlock, smemSize>>>(  p_P, p_O,p_dP, p_dO, p_chi, lambda, 
                                                    p_AA, p_Atb, Npix[0], Npix[1], Npix[2], blocks, single_probe);
            break;
        case 16:
           get_optimal_LSQ_step_ker< 16><<< dimGrid, dimBlock, smemSize>>>(  p_P, p_O,p_dP, p_dO, p_chi, lambda, 
                                                    p_AA, p_Atb, Npix[0], Npix[1], Npix[2], blocks, single_probe);
            break;
        case 8:
           get_optimal_LSQ_step_ker< 8><<< dimGrid, dimBlock, smemSize>>>(  p_P, p_O,p_dP, p_dO, p_chi, lambda, 
                                                    p_AA, p_Atb, Npix[0], Npix[1], Npix[2], blocks, single_probe);
            break;
        case 4:
           get_optimal_LSQ_step_ker< 4><<< dimGrid, dimBlock, smemSize>>>(  p_P, p_O,p_dP, p_dO, p_chi, lambda, 
                                                    p_AA, p_Atb, Npix[0], Npix[1], Npix[2], blocks, single_probe);
            break;
        case 2:
           get_optimal_LSQ_step_ker< 2><<< dimGrid, dimBlock, smemSize>>>(  p_P, p_O,p_dP, p_dO, p_chi, lambda, 
                                                    p_AA, p_Atb, Npix[0], Npix[1], Npix[2], blocks, single_probe);
            break;
        case 1:
           get_optimal_LSQ_step_ker< 1><<< dimGrid, dimBlock, smemSize>>>(  p_P, p_O,p_dP, p_dO, p_chi, lambda, 
                                                    p_AA, p_Atb, Npix[0], Npix[1], Npix[2], blocks, single_probe);
            break;
    }
  

    checkLastError("after kernel");           


    cudaThreadSynchronize();


    checkLastError("after kernel");

	// Wrap the result up as a MATLAB gpuArray for return. 
    plhs[0] = mxGPUCreateMxArrayOnGPU(m_AA);
	plhs[1] = mxGPUCreateMxArrayOnGPU(m_Atb);


    mxGPUDestroyGPUArray(m_P);
    mxGPUDestroyGPUArray(m_O);
    mxGPUDestroyGPUArray(m_dP);
    mxGPUDestroyGPUArray(m_dO);
    mxGPUDestroyGPUArray(m_chi);
    mxGPUDestroyGPUArray(m_AA);
    mxGPUDestroyGPUArray(m_Atb);




    
  return;
}




