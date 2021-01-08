
#include <algorithm>
#include <cuda_runtime_api.h>
#include "TV_texture.hpp"
#include <cuda.h>
#include <iostream>
#include <list>

#include "mex.h"
#include "gpu/mxGPUArray.h"

#define MAX(x,y) (x>y?x:y);
#define MIN(x,y) (x<y?x:y);
#define ABS(x) (x>0?x:-x);
#define INF (1023);

// limit of the number of neighbours 
static const unsigned MAX_NBRS = 256;
__constant__ float gC_neighbours_weights[MAX_NBRS];
__constant__ uint8_T gC_neighbours[MAX_NBRS];

typedef const unsigned int cuint; 
typedef const int cint; 

typedef texture<float, 3, cudaReadModeElementType> texture3D;

static texture3D ImgTexture, Xi_tex_0, Xi_tex_1, Xi_tex_2;

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

static bool bindVolumeDataTexture(const cudaArray* array, texture3D & Texture)
{
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	Texture.addressMode[0] = cudaAddressModeClamp;
	Texture.addressMode[1] = cudaAddressModeClamp;
	Texture.addressMode[2] = cudaAddressModeClamp;
	Texture.filterMode =  cudaFilterModeLinear; //cudaFilterModePoint
	Texture.normalized = false;

	cudaError err = cudaBindTextureToArray(Texture, array, channelDesc);
	checkLastError("cudaBindTextureToArray ");
	return true;
}

bool transferVolumeToArray(const float * img, cudaArray* array,  cuint X, cuint Y, cuint Z)
{

	//if (X * sizeof(float) > 2048) {
    //        mexPrintf("Volume is too large to be transfered to GPU array");
    //        return false;
	//}
	
    /* make volume array (no copying) */
     cudaPitchedPtr volume;
      volume.ptr = (float *)img; 
      volume.pitch = X * sizeof(float);
      volume.xsize = X;
      volume.ysize = Y;
     

        cudaExtent extent;
        extent.width = X;
        extent.height = Y;
        extent.depth = Z;

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
        sprintf(err, "total variation failed \n %s: %s. \n", msg, cudaGetErrorString(cudaStatus));
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
* SIMPLE GRADIENT BASED LOCAL TV
**/  


__global__ void kernel_local_TV3D_grad(float * p, const float tau, const float eps,
         cuint N, cuint M, cuint O, 
         cuint Xstart, cuint Ystart, cuint Zstart)  {

    // Location in a 3D matrix
    int m = Xstart+ blockIdx.x * blockDim.x + threadIdx.x;
    int n = Ystart+ blockIdx.y * blockDim.y + threadIdx.y;
    int o = Zstart+ blockIdx.z * blockDim.z + threadIdx.z;
 
    if (m < M & n < N & o < O)
    {
        int ms, ns,os, i,j,k, MN,  id ;
        cuint klm = 8; 
        MN = M*N;
        float tvg, nrm, p_mno; 
        float dp[3*klm];

        /********** gdv = div(grad( x );  ************/  

        for(i=0; i<2; i++) {
            ms = m+i;
            for(j=0; j<2; j++) {
                ns = n+j;
                for(k=0; k<2;k++)  {
                    id = 2*i+j+4*k;
                    // not all combinations are needed !!! 
                   if(id == 0 || id == 1 || id == 2 || id == 4) 
                    {
                    os = o+k;
                    p_mno =  tex3D(ImgTexture,ns+0.5f,ms+0.5f,os+0.5f);
                    dp[id] = p_mno - tex3D(ImgTexture,ns-1+0.5f,ms+0.5f,os+0.5f);
                    dp[id+klm] = p_mno - tex3D(ImgTexture,ns+0.5f,ms-1+0.5f,os+0.5f);
                    dp[id+2*klm] = p_mno - tex3D(ImgTexture,ns+0.5f,ms+0.5f,os-1+0.5f);
                    }
                }
            }
        }



        // get divergence of the gradient 
        tvg  = (dp[0] -  dp[1]) + (dp[klm] -  dp[2+klm])+ (dp[2*klm] -  dp[4+2*klm]);
        // get norm
        nrm = sqrt(dp[1]*dp[1] + dp[2+klm]*dp[2+klm] +  dp[4+2*klm]* dp[4+2*klm]);
        // avoid dividing by zero 
        nrm += eps;  

        //  update the image 
        p[m*N+n+o*MN] -= tau * tvg / nrm;
    }
}


/**
* CHAMBOLLE LOCAL TV SOLVER
**/



__global__ void kernel_local_TV3D_update_tex(float * p, 
        const float lambda, cuint N, cuint M, cuint O, 
         cuint Xstart, cuint Ystart, cuint Zstart)  {

    // Location in a 3D matrix
    int m = Xstart+ blockIdx.x * blockDim.x + threadIdx.x;
    int n = Ystart+ blockIdx.y * blockDim.y + threadIdx.y;
    int o = Zstart+ blockIdx.z * blockDim.z + threadIdx.z;

    if (m < M & n < N & o < O)
    {
        float div_xi;
        int ms, ns, os, MN, MNO;
        /***********  x = x - lambda*div( xi   ); % image update  ***********/
        ms = MAX(0,m-1);
        ns = MAX(0,n-1);
        os = MAX(0,o-1);
        MN = M*N; 
        MNO = M*N*O;
      // 10% of full iteration 
       div_xi  = tex3D(Xi_tex_0,n+0.5f,m+0.5f,o+0.5f) -  tex3D(Xi_tex_0,ns+0.5f,m+0.5f,o+0.5f) +
                 tex3D(Xi_tex_1,n+0.5f,m+0.5f,o+0.5f) -  tex3D(Xi_tex_1,n+0.5f,ms+0.5f,o+0.5f) +
                 tex3D(Xi_tex_2,n+0.5f,m+0.5f,o+0.5f) -  tex3D(Xi_tex_2,n+0.5f,m+0.5f,os+0.5f);
                 
    
      // 10% of the 3D TV total step 
       p[m*N+n+o*MN] -=  lambda*div_xi;   // save the updated reconstruction 
                                           // dont use the update this iteration 
    }
}

__global__ void kernel_local_TV3D_xi_tex(float * xi0,float * xi1,float * xi2,
        const float tau, const float lambda, cuint N, cuint M, cuint O, 
        cuint Xstart, cuint Ystart, cuint Zstart)  {

    // Location in a 3D matrix
    int m = Xstart+ blockIdx.x * blockDim.x + threadIdx.x;
    int n = Ystart+ blockIdx.y * blockDim.y + threadIdx.y;
    int o = Zstart+ blockIdx.z * blockDim.z + threadIdx.z;

    if (m < M & n < N & o < O)
    {
        int ms, ns,mss, nss,os, oss,  i,j,k, MN, MNO, id ; 
        MN = M*N;
        MNO = M*N*O;
        float d, p_mno = tex3D(ImgTexture,n+0.5f,m+0.5f,o+0.5f);
        float  gdv[3], div_xi8[8];


        // get new XI guess 
        /********** gdv = grad( div(xi) - x/lambda );  ************/  

        for(i=0; i<2; i++)
            {
            ms = MIN(m+i,M-1);
            mss = MAX(0,ms-1);
            for(j=0; j<2; j++) {
                ns = MIN(n+j,N-1);
                nss = MAX(0,ns-1);
                for(k=0; k<2;k++)  {
                    id = 2*i+j+4*k;
                    // not all cobinations are needed !!! 
                    if(id == 0 || id == 1 || id == 2 || id == 4) {
                    os = MIN(o+k,O-1);
                    oss = MAX(0,os-1);
                    div_xi8[2*i+j+4*k] = tex3D(Xi_tex_0,ns+0.5f,ms+0.5f,os+0.5f)  
                                        - tex3D(Xi_tex_0,nss+0.5f,ms+0.5f,os+0.5f) 
                                        + tex3D(Xi_tex_1,ns+0.5f,ms+0.5f,os+0.5f)  
                                        - tex3D(Xi_tex_1,ns+0.5f,mss+0.5f,os+0.5f)  
                                        + tex3D(Xi_tex_2,ns+0.5f,ms+0.5f,os+0.5f)  
                                        - tex3D(Xi_tex_2,ns+0.5f,mss+0.5f,oss+0.5f)  ;
                    }
                }
                }
            }

        ms = MIN(m+1,M-1);
        ns = MIN(n+1,N-1);
        os = MIN(o+1,O-1);
        gdv[0] = div_xi8[2]  -div_xi8[0]- (tex3D(ImgTexture,ns+0.5f,m+0.5f,o+0.5f)-p_mno)/lambda;
        gdv[1] = div_xi8[1]  -div_xi8[0]- (tex3D(ImgTexture,n+0.5f,ms+0.5f,o+0.5f)-p_mno)/lambda;
        gdv[2] = div_xi8[4]  -div_xi8[0]- (tex3D(ImgTexture,n+0.5f,m+0.5f,os+0.5f)-p_mno)/lambda;

       // for(i=0; i<3;i++)
         //  gdv_all[m*N+n+MN*o+i*MNO] = gdv[i];

        /**********
            %% isotropic 
            %d = sqrt(sum(gdv.^2));
            %% anisotropic 
            d = sum(abs(gdv));
        ****************/ 
        d = ABS(gdv[0]) 
        d += ABS(gdv[1]);
        d += ABS(gdv[2]);
        d = 1/( 1+tau*d );  

        /*****  chambolle step   **** */ 
        
        // around 50% of time needed to write the data back 
        xi0[m*N+n+o*MN] =  ( tex3D(Xi_tex_0,n+0.5f,m+0.5f,o+0.5f) + tau*gdv[0] )*d ; 
        xi1[m*N+n+o*MN] =  ( tex3D(Xi_tex_1,n+0.5f,m+0.5f,o+0.5f) + tau*gdv[1] )*d ; 
        xi2[m*N+n+o*MN] =  ( tex3D(Xi_tex_2,n+0.5f,m+0.5f,o+0.5f) + tau*gdv[2] )*d ; 

    }
}



/**
* NONLOCAL TV WEIGHTS 
**/ 

/*****************************
*******  TEXTURE BASED VERSIONS 
************************/

__global__ void kernel_nonlocal_weight_TV3D_tex( uint8_T *c, const uint8_T * neighbours, 
   cuint Nclose, cuint Nclose_min,  cuint Nnbrs, cint Rwin,  cint Rpatch, 
   cuint N, cuint M, cuint O, const float threshold, cuint Xstart, cuint Ystart, cuint Zstart)  {

    // Location in a 3D matrix
    int m = Xstart+ blockIdx.x * blockDim.x + threadIdx.x;
    int n = Ystart+ blockIdx.y * blockDim.y + threadIdx.y;
    int o = Zstart+ blockIdx.z * blockDim.z + threadIdx.z;

    if (m < M & n < N & o < O)
    {
        
    int const Nwin=Rwin*2+1;
    int const Nwin3 = Nwin*Nwin*Nwin;
    int const Nwin2 = Nwin*Nwin;



    int i,j,k,l, idx, idy, idz, idxx, idyy, idzz, ind_ijk;
    float D_sum,patch_diff,min_dist;

    float D_win[MAX_NBRS];
    for(i=0; i < MAX_NBRS; i++ )
        D_win[i] = INF; 
     
    int idxs, idys, idzs;
    for(i=-Rwin; i<=Rwin; i++) {
        for(j=-Rwin; j<=Rwin; j++)  {
            for (k=-Rwin; k<=Rwin; k++)  {
            D_sum = 0;
            for(idx = -Rpatch; idx<=Rpatch; idx++) 
                for(idy = -Rpatch; idy<=Rpatch; idy++)
                    for(idz = -Rpatch; idz<=Rpatch; idz++)
                    {

                    // original patch index 
                    idxx = m+idx; //idxx = MAX(0,idxx); idxx = MIN(M-1,idxx);
                    idyy = n+idy; //idyy = MAX(0,idyy); idyy = MIN(N-1,idyy);
                    idzz = o+idz; //idzz = MAX(0,idzz); idzz = MIN(O-1,idzz);
                    
                    idxs = idxx+i; //idxs = MAX(0,idxs); idxx = MIN(M-1,idxs); 
                    idys = idyy+j; //idys = MAX(0,idys); idyy = MIN(N-1,idys);
                    idzs = idzz+k; //idzz = MAX(0,idzs); idxx = MIN(O-1,idzs);
                    

                    // get difference between the shifted patches
                    // 70% of total time 
                    patch_diff = tex3D(ImgTexture,     idys+0.5f,  idxs+0.5f,idzs+0.5f) -
                                     tex3D(ImgTexture, idyy+  0.5f,  idxx+0.5f,  idzz+  0.5f);
                    // sum it up (ie convolution of squared data )
                     D_sum += patch_diff*patch_diff;

                    }
            ind_ijk = (k+Rwin)*Nwin2 + (i+Rwin)*Nwin + (j+Rwin);
            D_win[ind_ijk] = D_sum;
            }
        }
    }


    // remove the geometrically closest from consideration 
    for(k=0;k< Nnbrs ;k++)
        D_win[neighbours[k]-1] = INF;


    int min_ind=0;
    // find "most_similar" patches 
    for(k=0; k< Nclose; k++) {
        min_dist = INF;
        for(l=0; l<Nwin3; l++) { // find the closest patch
            min_dist = MIN(min_dist, D_win[l]);
            min_ind = (min_dist == D_win[l] ? l : min_ind);
        }
        D_win[min_ind] = INF;
        // save the position of the closest patch 
        if (min_dist  < threshold || k < Nclose_min)
            c[N*m + n + o*M*N + k*M*N*O] =  min_ind + 1;  // 20% of total time 
    }



    }
}



__global__ void kernel_nonlocal_TV3D_tex(float * p, const float *c_p0,
     const uint8_T *c,const uint8_T * neighbours, const float * Nbrs_weights,  const float dt, 
     const float eps, const float lambda, cuint Nclose,cuint Nnbrs,
     cint Rwin, cuint N, cuint M, cuint O, cuint Xstart, cuint Ystart, cuint Zstart)  {

    // Location in a 3D matrix
    int m = Xstart+ blockIdx.x * blockDim.x + threadIdx.x;
    int n = Ystart+ blockIdx.y * blockDim.y + threadIdx.y;
    int o = Zstart+ blockIdx.z * blockDim.z + threadIdx.z;


    if (m < M & n < N & o < O)
    {

        float dP = 0;
        int Nwin = 2 * Rwin + 1;
        int Nwin2 = Nwin*Nwin;
        int MNO = M*N*O;
        int ind_mno = n + N*m + M*N*o;
        float c_p_nmo = tex3D(ImgTexture, n+0.5f,m+0.5f,o+0.5f);
		float dp_tmp, adp_tmp;
		int i,j,k,l;
        int ind_close[MAX_NBRS];
        float weight;

        for (l=0; l < Nnbrs; l++) {
              ind_close[l] = gC_neighbours[l]-1;
              //weight[l] = (l == 0)? 1 : (0.1f);  // apply at least small weights to nbrs 
        }
        for (l=0; l < Nclose; l++) {
              ind_close[l+Nnbrs] = c[ind_mno + l*MNO]-1;
              //weight[l+Nnbrs] = gC_neighbours_weights[l];
        }

        for (l=0; l < Nnbrs  + Nclose ; l++)
        {   
            j = ind_close[l]%Nwin - Rwin;
            i = (ind_close[l]/Nwin)%Nwin - Rwin;
            k = (ind_close[l]/Nwin2) - Rwin;
            weight = (l == 0)? 1 : (0.1f);  // apply at least small weights to nbrs 
            weight = (l < Nnbrs ? weight : gC_neighbours_weights[l-Nnbrs]);  // weight differently in dependence on distance 
            weight = weight * (ind_close[l] >= 0); // ignore the too far nbrs 
            
            if (weight > 0) {  // read only if needed 
                dp_tmp = tex3D(ImgTexture, n+j+0.5f,m+i+0.5f,o+k+0.5f) - c_p_nmo;
                adp_tmp = ABS(dp_tmp);
                dP += weight * dp_tmp / (adp_tmp + eps); //  random noise dumping
            }
        }
        if (lambda > 0)
            p[ind_mno] += dt *( 2*dP + lambda*(c_p0[ind_mno] - c_p_nmo));
        else
            p[ind_mno] += dt * 2*dP; 


    }
}




/*****************
*******  TEXTURE BASED VERSIONS   2D
************************/



__global__ void kernel_nonlocal_weight_TV_tex( uint8_T *c, const uint8_T * neighbours, 
   cuint Nclose,  cuint Nnbrs, cint Rwin,  cint Rpatch, 
   cuint N, cuint M, const float threshold)  {

    // Location in a 2D matrix
    int m = blockIdx.x * blockDim.x + threadIdx.x;
    int n = blockIdx.y * blockDim.y + threadIdx.y;
    if (m < M & n < N)
    {
        
    int const Nwin=Rwin*2+1;
    int const Nwin2 = Nwin*Nwin;

    int i,j,k,l, idx, idy, idxx, idyy, Rwin_x, Rwin_y, Rpatch_x, Rpatch_y;
    float D_sum,patch_diff,min_dist;

    Rwin_x = MIN(Rwin, M-1);
    Rwin_y = MIN(Rwin, N-1);
    Rpatch_x = MIN(Rpatch, M-1);
    Rpatch_y = MIN(Rpatch, N-1);

    float D_win[MAX_NBRS];
    for(i=0; i < MAX_NBRS; i++ )
        D_win[i] = INF; 
   

    for(i=-Rwin_x; i<=Rwin_x; i++) {
        for(j=-Rwin_y; j<=Rwin_y; j++)  {
            D_sum = 0;
            for(idx = -Rpatch_x; idx<=Rpatch_x; idx++) {
                for(idy = -Rpatch_x; idy<=Rpatch_y; idy++)
                    {
                    // original patch index 
                    idxx = m+idx;
                    idyy = n+idy;
                    // get difference between the shifted patches
                    // presumably the slowest line 
                    patch_diff = tex3D(ImgTexture,     idyy+j+0.5f,idxx+i+0.5f,0.5f) -
                                     tex3D(ImgTexture, idyy+  0.5f,  idxx+0.5f,0.5f);
                    // sum it up (ie convolution of squared data )
                    D_sum += patch_diff*patch_diff;
                    }
                }
                D_win[(i+Rwin)*Nwin + j+Rwin] = D_sum;
            }
    }


    // remove the geometrically closest from consideration 
    for(k=0;k< Nnbrs ;k++)
        D_win[neighbours[k]-1] = INF;

    // find "most_similar" patches 
    for(k=0; k< Nclose; k++) {
        min_dist = INF;
        for(l=0; l<Nwin2; l++)  // find the closest patch
            min_dist = MIN(min_dist, D_win[l])
        for(l=0; l<Nwin2; l++)
            if(min_dist==D_win[l]) {
                D_win[l] = INF;
                break;
            }
        // save the position of the closest patch 
        c[N*m + n + k*M*N] = (min_dist  < threshold) ? l + 1 : 0; 
    }
    

    }
}



__global__ void kernel_nonlocal_TV_tex(float * p, const float * c_p, const float *c_p0,
     const uint8_T *c,const uint8_T * neighbours, const float dt, 
     const float eps, const float lambda, cuint Nclose,cuint Nnbrs,
     cint Rwin, cuint N, cuint M)  {

    // Location in a 2D matrix
    int m = blockIdx.x * blockDim.x + threadIdx.x;
    int n = blockIdx.y * blockDim.y + threadIdx.y;
    if (m < M & n < N)
    {

        float dP = 0;
        int Nwin = 2 * Rwin + 1;
        int ind_mn = n + N*m;
        float c_p_nm = tex3D(ImgTexture, n+0.5f,m+0.5f,0.5f);

		float dp_tmp, adp_tmp;
        int i,j,k;

        int ind_close[MAX_NBRS];
        for (k=0; k < Nnbrs; k++)
              ind_close[k] = neighbours[k]-1;
        for (k=0; k < Nclose; k++)
              ind_close[k+Nnbrs] = c[ind_mn + k*M*N]-1;

        for (k= 0 ; k < Nnbrs + Nclose ; k++)
        {   
            if (ind_close[k] == -1)
                continue;
            j = ind_close[k]%Nwin - Rwin;
            i = ind_close[k]/Nwin - Rwin;

            dp_tmp = tex3D(ImgTexture, n+j+0.5f,m+i+0.5f,0.5f) - c_p_nm;
			adp_tmp = ABS(dp_tmp);
			dP +=  dp_tmp / (adp_tmp + eps); //  random noise dumping
        }

        p[ind_mn] +=  dt *( 2*dP + lambda*(c_p0[ind_mn] - c_p_nm));
       // p[ind_mn] = dP; 
    }
}










/**
 * Host function called by MEX gateway. 
 */


    
/**
 * Nonlocal weights 
*/


void nonlocal_weight_TV_init( uint8_T *CloseInd, const float * Img,
       const uint8_T * neighbours,  cuint Nclose, cuint Nclose_min, cuint Nnbrs, cint Rwin,  cint Rpatch, 
       cuint N, cuint M, cuint O, const float threshold)
    {

    if (Nclose + Nnbrs > MAX_NBRS) {
        mexPrintf("Nclose + Nnbrs exceeded MAX_NBRS\n");
        return;
    }
    int Nwin = 2*Rwin+1;
    if ( (Nwin*Nwin*Nwin > MAX_NBRS & O > 1) |  (Nwin*Nwin > MAX_NBRS & O == 1) ) {
        mexPrintf("Nwin^3 exceeded MAX_NBRS\n");
        return;
    }
    if (M*N*O*4 > 1024e6) {  
        mexPrintf("Image size exceeded 1024MB, textures in  nonlocal TV will fail\n");
        return;
    }

    cudaMemcpyToSymbol(gC_neighbours, neighbours, Nnbrs*sizeof(uint8_T), 0, cudaMemcpyHostToDevice);

    /*  move image to the texture array  */
	cudaArray* cuArray = allocateVolumeArray(M,N,O);
    if (cuArray == 0)
        return;

	checkLastError("after allocateVolumeArray");
    transferVolumeToArray(Img, cuArray, M,N,O);
    checkLastError("after transferVolumeToArray\n \n ");
	bindVolumeDataTexture(cuArray, ImgTexture);




    if (O == 1) {
      // *************** 2-dim case ***************
      // Choose a reasonably sized number of threads in each dimension for the block.
        int const threadsPerBlockEachDim =  16;  // MAX THREAD is 1024 ~ 10*10*10 for 3D 
        dim3 const dimThread(threadsPerBlockEachDim, threadsPerBlockEachDim, 1);
        // Compute the thread block and grid sizes based on the board dimensions.
        int const blocksPerGrid_M = (M + threadsPerBlockEachDim - 1) / threadsPerBlockEachDim;
        int const blocksPerGrid_N = (N + threadsPerBlockEachDim - 1) / threadsPerBlockEachDim;
        dim3 dimBlock(blocksPerGrid_M, blocksPerGrid_N, 1);

         kernel_nonlocal_weight_TV_tex<<<dimBlock, dimThread>>>(CloseInd,neighbours,Nclose, Nnbrs, Rwin, Rpatch, M,N, threshold); 

    } else {
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

       // kernel_nonlocal_weight_TV3D_tex<<<dimBlock, dimThread>>>
       //         (CloseInd,neighbours,Nclose,Nclose_min, Nnbrs, Rwin, Rpatch, M,N,O, 
       //             threshold,0,0,0);

        std::list<cudaStream_t> streams;

        for ( int blockXstart=0; blockXstart < M; blockXstart += g_blockX)
            for ( int blockYstart=0; blockYstart < N; blockYstart += g_blockY)
                for ( int blockZstart=0; blockZstart < O; blockZstart += g_blockZ)
                    {
                    cudaStream_t stream;
                    cudaStreamCreate(&stream);
                    streams.push_back(stream);     
                    kernel_nonlocal_weight_TV3D_tex<<<dimBlock, dimThread, 0, stream>>>
                            (CloseInd,neighbours,Nclose,Nclose_min, Nnbrs, Rwin, Rpatch, M,N,O, 
                                threshold,blockXstart,blockYstart,blockZstart);                
                    }

        checkLastError("after kernel");           
        cudaThreadSynchronize(); 
        for (std::list<cudaStream_t>::iterator iter = streams.begin(); iter != streams.end(); ++iter)
                cudaStreamDestroy(*iter);

        streams.clear();
 
    }

    cudaTextForceKernelsCompletion();
 
	cudaFreeArray(cuArray);
	checkLastError("after cudaFreeArray");
	cudaUnbindTexture(ImgTexture);
	checkLastError("cudaUnbindTexture");

}


/**
 * nonlocal total variation 
*/

void nonlocal_TV_init( float * Img, const float * Img0, const uint8_T *CloseInd, 
       const uint8_T * neighbours, const float * Nbrs_weights,  
       const float dt, const float eps, const float lambda,cuint Nclose, cuint Nnbrs, cint Rwin, 
       cuint N, cuint M, cuint O,  cuint Niter)
    {

    if (Nclose + Nnbrs > MAX_NBRS) {
        mexPrintf("Nclose + Nnbrs exceeded MAX_NBRS\n");
        return;
    }
    int Nwin = 2*Rwin+1;
    if ( (Nwin*Nwin*Nwin > MAX_NBRS & O > 1) |  (Nwin*Nwin > MAX_NBRS & O == 1) ) {
        mexPrintf("Nwin^3 exceeded MAX_NBRS\n");
        return;
    }
    if (M*N*O*4 > 1024e6) {  
        mexPrintf("Image size exceeded 1024MB, textures in  nonlocal TV will fail\n");
        return;
    }

    cudaMemcpyToSymbol(gC_neighbours, neighbours, Nnbrs*sizeof(uint8_T), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(gC_neighbours_weights, Nbrs_weights, Nclose*sizeof(float), 0, cudaMemcpyHostToDevice);

    /*  move image to the texture array  */
	cudaArray* cuArray = allocateVolumeArray(M,N,O);
    if (cuArray == 0)
        return;
	checkLastError("after allocateVolumeArray");


    if (O == 1) {
        // *************** 2-dim case ***************
        // Choose a reasonably sized number of threads in each dimension for the block.
        int const threadsPerBlockEachDim =  16;  // MAX THREAD is 1024 ~ 10*10*10 for 3D 
        dim3 const dimThread(threadsPerBlockEachDim, threadsPerBlockEachDim, 1);
        // Compute the thread block and grid sizes based on the board dimensions.
        int const blocksPerGrid_M = (M + threadsPerBlockEachDim - 1) / threadsPerBlockEachDim;
        int const blocksPerGrid_N = (N + threadsPerBlockEachDim - 1) / threadsPerBlockEachDim;
        dim3 dimBlock(blocksPerGrid_M, blocksPerGrid_N, 1);

          for(int i=0 ; i < Niter; i++) {
                // mexPrintf("Iter %i\n", i);
                if (!transferVolumeToArray(Img, cuArray, M,N,O))
                    return;
                checkLastError("after transferVolumeToArray\n \n ");
                bindVolumeDataTexture(cuArray, ImgTexture);
                kernel_nonlocal_TV_tex<<<dimBlock, dimThread>>>(Img, Img0, Img0, CloseInd, neighbours, 
                                dt, eps,lambda, Nclose, Nnbrs, Rwin, M,N);
            }

    } else {
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

                for(int i=0 ; i < Niter; i++) {
                    // mexPrintf("Iter %i\n", i);
                    //fflush(stdout);

                    if (!transferVolumeToArray(Img, cuArray, M,N,O))
                        return;
                    checkLastError("after transferVolumeToArray\n \n ");
                    bindVolumeDataTexture(cuArray, ImgTexture);


                    for ( int blockXstart=0; blockXstart < M; blockXstart += g_blockX)
                        for ( int blockYstart=0; blockYstart < N; blockYstart += g_blockY)
                            for ( int blockZstart=0; blockZstart < O; blockZstart += g_blockZ) {
                                cudaStream_t stream;
                                cudaStreamCreate(&stream);
                                streams.push_back(stream);     
                                kernel_nonlocal_TV3D_tex<<<dimBlock, dimThread, 0, stream>>>(Img, Img0, 
                                        CloseInd, neighbours, Nbrs_weights,
                                        dt, eps,lambda, Nclose, Nnbrs, Rwin, M,N,O, 
                                        blockXstart,blockYstart,blockZstart);             
                                }

                    checkLastError("after kernel");           
                    cudaThreadSynchronize(); 
                	checkLastError("after cudaThreadSynchronize");

                    for (std::list<cudaStream_t>::iterator iter = streams.begin(); iter != streams.end(); ++iter)
                            cudaStreamDestroy(*iter);

                    streams.clear();
                }

            }

    cudaTextForceKernelsCompletion();
	cudaUnbindTexture(ImgTexture);
	checkLastError("cudaUnbindTexture");
	cudaFreeArray(cuArray);
	checkLastError("after cudaFreeArray");

}

/**
 * local total variation 
*/

void local_TV_init( float * Img, const float dt, const float eps, cuint M, cuint N, cuint O, cuint  Niter)
{
     if (M*N*O*4 > 1024e6) {  
        mexPrintf("Image size exceeded 1024MB, textures in  nonlocal TV will fail\n");
        return;
    }
    /*  move image to the texture array  */
	cudaArray* cuArray = allocateVolumeArray(M,N,O);
    if (cuArray == 0)
        return;
	checkLastError("after allocateVolumeArray");

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

    for(int i=0 ; i < Niter; i++) {
        if (!transferVolumeToArray(Img, cuArray, M,N,O))
            return;
        checkLastError("after transferVolumeToArray\n \n ");
        bindVolumeDataTexture(cuArray, ImgTexture);

        for ( int blockXstart=0; blockXstart < M; blockXstart += g_blockX)
            for ( int blockYstart=0; blockYstart < N; blockYstart += g_blockY)
                for ( int blockZstart=0; blockZstart < O; blockZstart += g_blockZ) {
                    cudaStream_t stream;
                    cudaStreamCreate(&stream);
                    streams.push_back(stream);     
                    kernel_local_TV3D_grad<<<dimBlock, dimThread>>>(Img, dt, eps, M,N,O,
                             blockXstart,blockYstart,blockZstart);             
                    }
        checkLastError("after kernel");           
        cudaThreadSynchronize(); 
        checkLastError("after cudaThreadSynchronize");
        for (std::list<cudaStream_t>::iterator iter = streams.begin(); iter != streams.end(); ++iter)
                cudaStreamDestroy(*iter);
        streams.clear();
    }


    cudaTextForceKernelsCompletion();
    cudaUnbindTexture(ImgTexture);
    checkLastError("cudaUnbindTexture");
    cudaFreeArray(cuArray);
    checkLastError("after cudaFreeArray");
}



void local_TV_chambolle_init( float * Img, float ** Xi, const float dt, const float tau, 
                                cuint M, cuint N, cuint O, cuint  Niter)
{
     if (M*N*O*4 > 1024e6) {  
        mexPrintf("Image size exceeded 1024MB, textures in  nonlocal TV will fail\n");
        return;
    }

    // Choose a reasonably sized number of threads in each dimension for the block.
    int const threadsPerBlockEachDim =  10;  // MAX THREAD is 1024 ~ 10*10*10 for 3D 
    // Compute the thread block and grid sizes based on the board dimensions.
    int const blocksPerGrid_M = (M + threadsPerBlockEachDim - 1) / threadsPerBlockEachDim;
    int const blocksPerGrid_N = (N + threadsPerBlockEachDim - 1) / threadsPerBlockEachDim;
    int const blocksPerGrid_O = (O + threadsPerBlockEachDim - 1) / threadsPerBlockEachDim;

    dim3 const dimBlock(blocksPerGrid_M, blocksPerGrid_N, blocksPerGrid_O);
    dim3 const dimThread(threadsPerBlockEachDim, threadsPerBlockEachDim, (O > 1 ? threadsPerBlockEachDim:1));

   // mexPrintf("Threads %i %i %i \n ", blocksPerGrid_M, blocksPerGrid_N, blocksPerGrid_O);
   // mexPrintf("Blocks %i %i %i \n ", dimThread.x, dimThread.y, dimThread.z);


   /*  move image to the texture array  */
	cudaArray * cuArrayImg = allocateVolumeArray(M,N,O);
    if (cuArrayImg == 0)
        return;
    
	cudaArray * cuArrayXi[3];
    for (int i = 0; i < 3; i++) {
        cuArrayXi[i] = allocateVolumeArray(M,N,O);
        if (cuArrayXi[i] == 0)
            return;
    }
	if (!checkLastError("after allocateVolumeArray"))
        return;

    for (int i = 0; i < 3; i++)
        transferVolumeToArray(Xi[i], cuArrayXi[i], M,N,O);

    bindVolumeDataTexture(cuArrayXi[0], Xi_tex_0);
    bindVolumeDataTexture(cuArrayXi[1], Xi_tex_1);
    bindVolumeDataTexture(cuArrayXi[2], Xi_tex_2);


    for(int iter=0 ; iter < Niter; iter++) {
        //mexPrintf("Iter %i\n", iter);
        // fflush(stdout);

        if (!transferVolumeToArray(Img, cuArrayImg, M,N,O))
            return;       
        bindVolumeDataTexture(cuArrayImg, ImgTexture);
        if (!checkLastError("after bindVolumeDataTexture"))
            return;

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



        if (!transferVolumeToArray(Img, cuArrayImg, M,N,O))
            return;
        checkLastError("after transferVolumeToArray\n \n ");
        bindVolumeDataTexture(cuArrayImg, ImgTexture);

        // ******** calculate the chambolle dual space Xi *************
        for ( int blockXstart=0; blockXstart < M; blockXstart += g_blockX)
            for ( int blockYstart=0; blockYstart < N; blockYstart += g_blockY)
                for ( int blockZstart=0; blockZstart < O; blockZstart += g_blockZ) {
                    cudaStream_t stream;
                    cudaStreamCreate(&stream);
                    streams.push_back(stream); 
                    kernel_local_TV3D_xi_tex<<<dimBlock, dimThread, 0, stream>>>(Xi[0],Xi[1],Xi[2],tau, dt,
                    N,M,O, blockXstart,blockYstart,blockZstart);
                    }

        checkLastError("after kernel");           

        for (std::list<cudaStream_t>::iterator iter = streams.begin(); iter != streams.end(); ++iter)
                cudaStreamDestroy(*iter);

        streams.clear();

        cudaTextForceKernelsCompletion();

        if (!checkLastError("after kernel_local_TV3D_xi_tex"))
            return;
        for (int i = 0; i < 3; i++)
            transferVolumeToArray(Xi[i], cuArrayXi[i], M,N,O);
        if (!checkLastError("after transferVolumeToArray"))
            return;

        bindVolumeDataTexture(cuArrayXi[0], Xi_tex_0);
        bindVolumeDataTexture(cuArrayXi[1], Xi_tex_1);
        bindVolumeDataTexture(cuArrayXi[2], Xi_tex_2);
      
        if (!checkLastError("after bindVolumeDataTexture"))
            return;

        // ******** calculate the chambolle real space update *************

         for ( int blockXstart=0; blockXstart < M; blockXstart += g_blockX)
                for ( int blockYstart=0; blockYstart < N; blockYstart += g_blockY)
                    for ( int blockZstart=0; blockZstart < O; blockZstart += g_blockZ) {
                        cudaStream_t stream;
                        cudaStreamCreate(&stream);
                        streams.push_back(stream); 
                        kernel_local_TV3D_update_tex<<<dimBlock, dimThread, 0, stream>>>
                            (Img, dt, M,N,O, blockXstart,blockYstart,blockZstart);
                        }

        checkLastError("after kernel");           
        cudaThreadSynchronize(); 
        for (std::list<cudaStream_t>::iterator iter = streams.begin(); iter != streams.end(); ++iter)
                cudaStreamDestroy(*iter);

        streams.clear();
        cudaTextForceKernelsCompletion();

    }

    cudaTextForceKernelsCompletion();
    cudaUnbindTexture(ImgTexture);
    cudaUnbindTexture(Xi_tex_0);
    cudaUnbindTexture(Xi_tex_1);
    cudaUnbindTexture(Xi_tex_2);


	checkLastError("after cudaThreadSynchronize");
	cudaFreeArray(cuArrayImg);
    for (int i = 0; i < 3; i++)
       cudaFreeArray(cuArrayXi[i]);

	checkLastError("after cudaFreeArray");
	checkLastError("cudaUnbindTexture");

}



    

