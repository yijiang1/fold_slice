
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




#include <algorithm>
#include <cuda_runtime_api.h>
#include "interp3_gpu.hpp"
#include <cuda.h>
#include <iostream>
#include <list>

#include "mex.h"
#include "gpu/mxGPUArray.h"

#define MAX(x,y) (x>y?x:y);
#define MIN(x,y) (x<y?x:y);
#define ABS(x) (x>0?x:-x);
#define INF (1023);

typedef const unsigned int cuint; 
typedef const int cint; 

typedef texture<float, 3, cudaReadModeElementType> texture3D;

static texture3D ImgTexture, X_tex, Y_tex, Z_tex;

// splitting volume on smaller blocks to prevent GPU crashes 
static cuint g_blockX = 256;
static cuint g_blockY = 256;
static cuint g_blockZ = 256;


cudaArray* allocateVolumeArray( cuint X, cuint Y, cuint Z)
{
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    cudaArray* cuArray;
    cudaExtent extent;
    extent.width = X;
    extent.height = Y;
    extent.depth = Z;
    cudaError err = cudaMalloc3DArray(&cuArray, &channelDesc, extent);
    if (err != cudaSuccess) {
        mexPrintf  ("Failed to allocate %dx%dx%d GPU array\n",X,Y,Z);
        return 0;
    }
    return cuArray;
}

static bool bindVolumeDataTexture(const cudaArray* array, texture3D & Texture, bool normalized)
{
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	Texture.addressMode[0] = cudaAddressModeClamp;
	Texture.addressMode[1] = cudaAddressModeClamp;
	Texture.addressMode[2] = cudaAddressModeClamp;
	Texture.filterMode =  cudaFilterModeLinear; //cudaFilterModePoint
	Texture.normalized = normalized;

	cudaError err = cudaBindTextureToArray(Texture, array, channelDesc);
	checkLastError("cudaBindTextureToArray ");
	return true;
}

bool transferVolumeToArray(const mxGPUArray * m_img, cudaArray *& array)
{
    
    mwSize  const * dimensions = mxGPUGetDimensions(m_img);
    mwSize Ndim =  mxGPUGetNumberOfDimensions(m_img);
    int M = (int)dimensions[0];
    int N = (int)dimensions[1];
    int O = Ndim > 2 ? (int)dimensions[2] :  1;
    
    // get the values into float array 
    const float * img =(const float *)mxGPUGetDataReadOnly(m_img);

    array = allocateVolumeArray(M,N,O);
    if (array == 0)
        return false;

	if (M * sizeof(float) > 2048) {
            mexPrintf("Volume is too large to be transfered to GPU array");
            return false;
	}
	
    /* make volume array (no copying) */
        cudaPitchedPtr volume;
        volume.ptr = (float *)img; 
        volume.pitch = M * sizeof(float);
        volume.xsize = M;
        volume.ysize = N;

        cudaExtent extent;
        extent.width = M;
        extent.height = N;
        extent.depth = O;

        cudaMemcpy3DParms p;
        cudaPos zp = { 0, 0, 0 };
        p.srcArray = 0;
        p.srcPos = zp;
        p.srcPtr = volume;
        p.dstArray = array;
        p.dstPtr.ptr = 0;
        p.dstPtr.pitch = 0;
        p.dstPtr.xsize = 0;
        p.dstPtr.ysize = 0;
        p.dstPos = zp;
        p.extent = extent;
        p.kind = cudaMemcpyDeviceToDevice;

        cudaError err = cudaMemcpy3D(&p);

        if (!checkLastError("transferVolumeToArray cudaMemcpy3D"))
            return false;
            
        return true;
}

int checkLastError(char * msg)
{
    cudaError_t cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        char err[512];
        sprintf(err, "interp3 variation failed \n %s: %s. \n", msg, cudaGetErrorString(cudaStatus));
        mexErrMsgTxt(err);
        return 0;
    }
    return 1;
}
bool cudaTextForceKernelsCompletion()
{
        cudaError_t returnedCudaError = cudaThreadSynchronize();
        if (returnedCudaError != cudaSuccess) {
                fprintf(stderr, "Failed to force completion of cuda kernels: %d: %s. \n ", returnedCudaError, cudaGetErrorString(returnedCudaError));
                return false;
        }
        return true;
}

/** 
* TEXTURE TRILINEAR INTERPOLATION 
**/  


__global__ void kernel_interp3(float * p, cuint N, cuint M, cuint O, 
         cuint Xstart, cuint Ystart, cuint Zstart)  {

    // Location in a 3D matrix
    mwSize m = Xstart+ blockIdx.x * blockDim.x + threadIdx.x;
    mwSize n = Ystart+ blockIdx.y * blockDim.y + threadIdx.y;
    mwSize o = Zstart+ blockIdx.z * blockDim.z + threadIdx.z;
 
    if (m < M & n < N & o < O)
    {
        float xs, ys, zs;  // shifted coordinates 
        float  mn, nn, on; // normalized coordinates 
        mn = (float)(m)/M; 
        nn = (float)(n)/N; 
        on = (float)(o)/O; 
      //  mn = (m+0.5f); 
      //  nn = (n+0.5f); 
      //  on = (o+0.5f);
        // load deformed coordinates 
        xs =  m+0.5f - tex3D(X_tex,mn, nn, on);
        ys =  n+0.5f - tex3D(Y_tex,mn, nn, on);
        zs =  o+0.5f - tex3D(Z_tex,mn, nn, on);
        // get trilinear interplation
        bool outsiders = (xs > 0) & (ys > 0) & (zs > 0) &
                            (xs < M) & (ys < N) & (zs < O);
       float p_val = (outsiders ? tex3D(ImgTexture,xs,ys,zs) : 0);
        //  write the interpolation to the output 
        p[(n)*N+(m)+(o)*M*N] = p_val;
    }
}


/**
 * Host function called by MEX gateway. 
 */


   

void  interp3_init( float * p, const mxGPUArray * p0, const mxGPUArray *m_X, const mxGPUArray *m_Y, const mxGPUArray *m_Z, cuint M, cuint N, cuint O)
{
    if (M*N*O*4 > 1024e6) {  
        mexPrintf("Image size exceeded 1024MB, textures in interp3 will fail\n");
        return;
    }

    /*  move image to the texture array  */
	cudaArray* cuArray, *cuArrayX, *cuArrayY, *cuArrayZ; 
	checkLastError("after allocateVolumeArray");
    transferVolumeToArray(p0, cuArray);
    checkLastError("after transferVolumeToArray\n \n ");
	bindVolumeDataTexture(cuArray, ImgTexture, false);

    /*  move X deformation to the texture array  */
	checkLastError("after allocateVolumeArray");
    transferVolumeToArray(m_X, cuArrayX);
    checkLastError("after transferVolumeToArray\n \n ");
	bindVolumeDataTexture(cuArrayX, X_tex, true);

    /*  move Y deformation to the texture array  */
	checkLastError("after allocateVolumeArray");
    transferVolumeToArray(m_Y, cuArrayY);
    checkLastError("after transferVolumeToArray\n \n ");
	bindVolumeDataTexture(cuArrayY, Y_tex, true);

     /*  move Z deformation to the texture array  */
	checkLastError("after allocateVolumeArray");
    transferVolumeToArray(m_Z, cuArrayZ);
    checkLastError("after transferVolumeToArray\n \n ");
	bindVolumeDataTexture(cuArrayZ, Z_tex, true);






     // *************** 3-dim case ***************
    // Choose a reasonably sized number of threads in each dimension for the block.
    int const threadsPerBlockEachDim =  10;  // MAX THREAD is 1024 ~ 10*10*10 for 3D 
    dim3 const dimThread(threadsPerBlockEachDim, threadsPerBlockEachDim, threadsPerBlockEachDim);
    //mexPrintf("Thread %i %i %i \n ", dimThread.x, dimThread.y, dimThread.z);
    // Compute the thread block and grid sizes based on the board dimensions.
    int const blocksPerGrid_M = (g_blockX + threadsPerBlockEachDim - 1) / threadsPerBlockEachDim;
    int const blocksPerGrid_N = (g_blockY + threadsPerBlockEachDim - 1) / threadsPerBlockEachDim;
    int const blocksPerGrid_O = (g_blockZ + threadsPerBlockEachDim - 1) / threadsPerBlockEachDim;
    dim3 dimBlock(blocksPerGrid_M, blocksPerGrid_N, blocksPerGrid_O);
    //mexPrintf("Block %i %i %i \n ", blocksPerGrid_M, blocksPerGrid_N, blocksPerGrid_O);

    std::list<cudaStream_t> streams;

    for ( int blockXstart=0; blockXstart < M; blockXstart += g_blockX)
        for ( int blockYstart=0; blockYstart < N; blockYstart += g_blockY)
            for ( int blockZstart=0; blockZstart < O; blockZstart += g_blockZ)
                {
                cudaStream_t stream;
                cudaStreamCreate(&stream);
                streams.push_back(stream);     
                kernel_interp3<<<dimBlock, dimThread, 0, stream>>>
                        (p,M,N,O, blockXstart,blockYstart,blockZstart);                
                }

    checkLastError("after kernel");           
    cudaThreadSynchronize(); 
    for (std::list<cudaStream_t>::iterator iter = streams.begin(); iter != streams.end(); ++iter)
            cudaStreamDestroy(*iter);

    streams.clear();
 


    // clear memory , unbind textures 

    cudaTextForceKernelsCompletion();

	cudaFreeArray(cuArray);
	cudaFreeArray(cuArrayX);
	cudaFreeArray(cuArrayY);
	cudaFreeArray(cuArrayZ);

	checkLastError("after cudaFreeArray");
	cudaUnbindTexture(ImgTexture);
	cudaUnbindTexture(X_tex);
	cudaUnbindTexture(Y_tex);
	cudaUnbindTexture(Z_tex);
	checkLastError("cudaUnbindTexture");



}



    

