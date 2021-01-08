/*
-----------------------------------------------------------------------
Copyright: 2010-2015, iMinds-Vision Lab, University of Antwerp
           2014-2015, CWI, Amsterdam

Contact: astra@uantwerpen.be
Website: http://sf.net/projects/astra-toolbox

This file is part of the ASTRA Toolbox.


The ASTRA Toolbox is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The ASTRA Toolbox is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the ASTRA Toolbox. If not, see <http://www.gnu.org/licenses/>.

-----------------------------------------------------------------------
$Id$
*/

#include <cstdio>
#include <cassert>
#include "util3d.h"
#include <ctime>

#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

//#include "../2d/util.h"

#include "astra/Logging.h"
#include "mex.h"


namespace astraCUDA3d {


	cudaPitchedPtr allocateVolumeData(const SDimensions3D& dims)
	{
		cudaExtent extentV;
		extentV.width = dims.iVolX*sizeof(float);
		extentV.height = dims.iVolY;
		extentV.depth = dims.iVolZ;

		cudaPitchedPtr volData;

		cudaError err = cudaMalloc3D(&volData, extentV);
		if (err != cudaSuccess) {
			astraCUDA3d::reportCudaError(err);
			ASTRA_ERROR("Failed to allocate %dx%dx%d GPU buffer", dims.iVolX, dims.iVolY, dims.iVolZ);
			volData.ptr = 0;
			// TODO: return 0 somehow?
		}

		return volData;
	}
	cudaPitchedPtr allocateProjectionData(const SDimensions3D& dims)
	{
		cudaExtent extentP;
		extentP.width = dims.iProjU*sizeof(float);
		extentP.height = dims.iProjAngles;
		extentP.depth = dims.iProjV;

		cudaPitchedPtr projData;

		cudaError err = cudaMalloc3D(&projData, extentP);
		if (err != cudaSuccess) {
			mexPrintf("Failed to allocate %dx%dx%d GPU buffer", dims.iProjU, dims.iProjAngles, dims.iProjV);
			projData.ptr = 0;
			// TODO: return 0 somehow?
		}

		return projData;
	}
	bool zeroVolumeData(cudaPitchedPtr& D_data, const SDimensions3D& dims)
	{
		char* t = (char*)D_data.ptr;
		cudaError err;

		for (unsigned int z = 0; z < dims.iVolZ; ++z) {
			err = cudaMemset2D(t, D_data.pitch, 0, dims.iVolX*sizeof(float), dims.iVolY);
			ASTRA_CUDA_ASSERT(err);
			t += D_data.pitch * dims.iVolY;
		}
		return true;
	}
	bool zeroProjectionData(cudaPitchedPtr& D_data, const SDimensions3D& dims)
	{
		char* t = (char*)D_data.ptr;
		cudaError err;

		for (unsigned int z = 0; z < dims.iProjV; ++z) {
			err = cudaMemset2D(t, D_data.pitch, 0, dims.iProjU*sizeof(float), dims.iProjAngles);
			ASTRA_CUDA_ASSERT(err);
			t += D_data.pitch * dims.iProjAngles;
		}

		return true;
	}
	bool copyVolumeToDevice(const float* data, cudaPitchedPtr& D_data, const SDimensions3D& dims, unsigned int pitch)
	{
		if (!pitch)
			pitch = dims.iVolX;

		cudaPitchedPtr ptr;
		ptr.ptr = (void*)data; // const cast away
		ptr.pitch = pitch*sizeof(float);
		ptr.xsize = dims.iVolX*sizeof(float);
		ptr.ysize = dims.iVolY;

		cudaExtent extentV;
		extentV.width = dims.iVolX*sizeof(float);
		extentV.height = dims.iVolY;
		extentV.depth = dims.iVolZ;

		cudaPos zp = { 0, 0, 0 };

		cudaMemcpy3DParms p;
		p.srcArray = 0;
		p.srcPos = zp;
		p.srcPtr = ptr;
		p.dstArray = 0;
		p.dstPos = zp;
		p.dstPtr = D_data;
		p.extent = extentV;
		p.kind = cudaMemcpyHostToDevice;

		cudaError err;
		err = cudaMemcpy3D(&p);
		ASTRA_CUDA_ASSERT(err);

		return err == cudaSuccess;
	}

	bool copyProjectionsToDevice(const float* data, cudaPitchedPtr& D_data, const SDimensions3D& dims, unsigned int pitch)
	{
		if (!pitch)
			pitch = dims.iProjU;

		cudaPitchedPtr ptr;
		ptr.ptr = (void*)data; // const cast away
		ptr.pitch = pitch*sizeof(float);
		ptr.xsize = dims.iProjU*sizeof(float);
		ptr.ysize = dims.iProjAngles;

		cudaExtent extentV;
		extentV.width = dims.iProjU*sizeof(float);
		extentV.height = dims.iProjAngles;
		extentV.depth = dims.iProjV;

		cudaPos zp = { 0, 0, 0 };

		cudaMemcpy3DParms p;
		p.srcArray = 0;
		p.srcPos = zp;
		p.srcPtr = ptr;
		p.dstArray = 0;
		p.dstPos = zp;
		p.dstPtr = D_data;
		p.extent = extentV;
		p.kind = cudaMemcpyHostToDevice;

		cudaError err;
		err = cudaMemcpy3D(&p);
		ASTRA_CUDA_ASSERT(err);

		return err == cudaSuccess;
	}

	bool copyVolumeFromDevice(float* data, const cudaPitchedPtr& D_data, const SDimensions3D& dims, unsigned int pitch)
	{
		if (!pitch)
			pitch = dims.iVolX;

		cudaPitchedPtr ptr;
		ptr.ptr = data;
		ptr.pitch = pitch*sizeof(float);
		ptr.xsize = dims.iVolX*sizeof(float);
		ptr.ysize = dims.iVolY;

		cudaExtent extentV;
		extentV.width = dims.iVolX*sizeof(float);
		extentV.height = dims.iVolY;
		extentV.depth = dims.iVolZ;

		cudaPos zp = { 0, 0, 0 };

		cudaMemcpy3DParms p;
		p.srcArray = 0;
		p.srcPos = zp;
		p.srcPtr = D_data;
		p.dstArray = 0;
		p.dstPos = zp;
		p.dstPtr = ptr;
		p.extent = extentV;
		p.kind = cudaMemcpyDeviceToHost;

		cudaError err;
		err = cudaMemcpy3D(&p);
		ASTRA_CUDA_ASSERT(err);

		return err == cudaSuccess;
	}
	bool copyProjectionsFromDevice(float* data, const cudaPitchedPtr& D_data, const SDimensions3D& dims, unsigned int pitch)
	{
		if (!pitch)
			pitch = dims.iProjU;

		cudaPitchedPtr ptr;
		ptr.ptr = data;
		ptr.pitch = pitch*sizeof(float);
		ptr.xsize = dims.iProjU*sizeof(float);
		ptr.ysize = dims.iProjAngles;

		cudaExtent extentV;
		extentV.width = dims.iProjU*sizeof(float);
		extentV.height = dims.iProjAngles;
		extentV.depth = dims.iProjV;

		cudaPos zp = { 0, 0, 0 };

		cudaMemcpy3DParms p;
		p.srcArray = 0;
		p.srcPos = zp;
		p.srcPtr = D_data;
		p.dstArray = 0;
		p.dstPos = zp;
		p.dstPtr = ptr;
		p.extent = extentV;
		p.kind = cudaMemcpyDeviceToHost;

		cudaError err;
		err = cudaMemcpy3D(&p);
		ASTRA_CUDA_ASSERT(err);

		return err == cudaSuccess;
	}

	bool duplicateVolumeData(cudaPitchedPtr& D_dst, const cudaPitchedPtr& D_src, const SDimensions3D& dims)
	{
		cudaExtent extentV;
		extentV.width = dims.iVolX*sizeof(float);
		extentV.height = dims.iVolY;
		extentV.depth = dims.iVolZ;

		cudaPos zp = { 0, 0, 0 };

		cudaMemcpy3DParms p;
		p.srcArray = 0;
		p.srcPos = zp;
		p.srcPtr = D_src;
		p.dstArray = 0;
		p.dstPos = zp;
		p.dstPtr = D_dst;
		p.extent = extentV;
		p.kind = cudaMemcpyDeviceToDevice;

		cudaError err;
		err = cudaMemcpy3D(&p);
		ASTRA_CUDA_ASSERT(err);

		return err == cudaSuccess;
	}
	bool duplicateProjectionData(cudaPitchedPtr& D_dst, const cudaPitchedPtr& D_src, const SDimensions3D& dims)
	{
		cudaExtent extentV;
		extentV.width = dims.iProjU*sizeof(float);
		extentV.height = dims.iProjAngles;
		extentV.depth = dims.iProjV;

		cudaPos zp = { 0, 0, 0 };

		cudaMemcpy3DParms p;
		p.srcArray = 0;
		p.srcPos = zp;
		p.srcPtr = D_src;
		p.dstArray = 0;
		p.dstPos = zp;
		p.dstPtr = D_dst;
		p.extent = extentV;
		p.kind = cudaMemcpyDeviceToDevice;

		cudaError err;
		err = cudaMemcpy3D(&p);
		ASTRA_CUDA_ASSERT(err);

		return err == cudaSuccess;
	}



	// TODO: Consider using a single array of size max(proj,volume) (per dim)
	//       instead of allocating a new one each time

	cudaArray* allocateVolumeArray(const SDimensions3D& dims)
	{
		cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
		cudaArray* cuArray;
		cudaExtent extentA;
		extentA.width = dims.iVolX;
		extentA.height = dims.iVolY;
		extentA.depth = dims.iVolZ;
		cudaError err = cudaMalloc3DArray(&cuArray, &channelDesc, extentA);
		if (err != cudaSuccess) {
			mexPrintf("Failed to allocate %dx%dx%d GPU array", dims.iVolX, dims.iVolY, dims.iVolZ);
			return 0;
		}

		return cuArray;
	}
	cudaArray* allocateProjectionArray(const SDimensions3D& dims)
	{
		cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
		cudaArray* cuArray;
		cudaExtent extentA;
		extentA.width = dims.iProjU;
		extentA.height = dims.iProjAngles;
		extentA.depth = dims.iProjV;
		cudaError err = cudaMalloc3DArray(&cuArray, &channelDesc, extentA);

		if (err != cudaSuccess) {
			mexPrintf("Failed to allocate %dx%dx%d GPU array", dims.iProjU, dims.iProjAngles, dims.iProjV);
			return 0;
		}

		return cuArray;
	}

    bool bindDataTexture(const cudaArray* array, texture3D & Texture, cudaTextureAddressMode bordermode, bool normalized)
    {
        cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
        Texture.addressMode[0] = bordermode;
        Texture.addressMode[1] = bordermode;
        Texture.addressMode[2] = bordermode;
        Texture.filterMode = cudaFilterModeLinear;
        Texture.normalized = normalized;

        cudaError err = cudaBindTextureToArray(Texture, array, channelDesc);

        checkLastError("cudaBindTextureToArray cudaMemcpy3D");
        ASTRA_CUDA_ASSERT(err);

        //mexPrintf("Max texture size !!!  %i %i %i", cudaDeviceProp.maxTexture3D[0], cudaDeviceProp.maxTexture3D[1], cudaDeviceProp.maxTexture3D[2]);
        return true;
    }

    cudaArray *  transferDeformationToArray(const mxGPUArray * m_img)
    {
        mwSize  const * dimensions = mxGPUGetDimensions(m_img);
        mwSize Ndim =  mxGPUGetNumberOfDimensions(m_img);
        int M = (int)dimensions[0];
        int N = (int)dimensions[1];
        int O = Ndim > 2 ? (int)dimensions[2] :  1;

        SDimensions3D dims;
        dims.iVolX = M;
        dims.iVolY = N;
        dims.iVolZ = O;

        //mexPrintf("Deformation field size: %i %i %i \n", M,N,O); 


    	cudaArray* array = allocateVolumeArray(dims);

        // get the values into float array 
        const float * img =(const float *)mxGPUGetDataReadOnly(m_img);

        if (array == 0)
            return 0;

        if (M * sizeof(float) > 2048) {
			mexPrintf("Volume is too large to be transfered to GPU array");
            return 0;
        }

        // make volume array (no copying) 
        cudaPitchedPtr volume;
        volume.ptr = (float *)img; 
        volume.pitch = M * sizeof(float);
        volume.xsize = M;
        volume.ysize = N;
        

        transferVolumeToArray(volume, array,dims);

      // if (!checkLastError("transferDeformToArray cudaMemcpy3D"))
      //     return false;

        return array;
    }


	bool transferVolumeToArray(cudaPitchedPtr D_volumeData, cudaArray* array, const SDimensions3D& dims)
	{
		cudaExtent extentA;
		extentA.width = dims.iVolX;
		extentA.height = dims.iVolY;
		extentA.depth = dims.iVolZ;

		cudaMemcpy3DParms p;
		cudaPos zp = { 0, 0, 0 };
		p.srcArray = 0;
		p.srcPos = zp;
		p.srcPtr = D_volumeData;
		p.dstArray = array;
		p.dstPtr.ptr = 0;
		p.dstPtr.pitch = 0;
		p.dstPtr.xsize = 0;
		p.dstPtr.ysize = 0;
		p.dstPos = zp;
		p.extent = extentA;
		p.kind = cudaMemcpyDeviceToDevice;

		cudaError err = cudaMemcpy3D(&p);

		checkLastError("transferVolumeToArray cudaMemcpy3D"); 
		ASTRA_CUDA_ASSERT(err);
		// TODO: check errors
		return true;
	}

	bool transferProjectionsToArray(cudaPitchedPtr D_projData, cudaArray* array, const SDimensions3D& dims)
	{
		cudaExtent extentA;
		extentA.width = dims.iProjU;
		extentA.height = dims.iProjAngles;
		extentA.depth = dims.iProjV;

		cudaMemcpy3DParms p;
		cudaPos zp = { 0, 0, 0 };
		p.srcArray = 0;
		p.srcPos = zp;
		p.srcPtr = D_projData;
		p.dstArray = array;
		p.dstPtr.ptr = 0;
		p.dstPtr.pitch = 0;
		p.dstPtr.xsize = 0;
		p.dstPtr.ysize = 0;
		p.dstPos = zp;
		p.extent = extentA;
		p.kind = cudaMemcpyDeviceToDevice;

		cudaError err = cudaMemcpy3D(&p);
		checkLastError("transferProjectionsToArray cudaMemcpy3D");

		ASTRA_CUDA_ASSERT(err);

		// TODO: check errors

		return true;
	}


	bool cudaTextForceKernelsCompletion()
	{
		cudaError_t returnedCudaError = cudaThreadSynchronize();

		if (returnedCudaError != cudaSuccess) {
			//FIXME 
			fprintf(stderr, "Failed to force completion of cuda kernels: %d: %s. \n ", returnedCudaError, cudaGetErrorString(returnedCudaError));
			ASTRA_ERROR("Failed to force completion of cuda kernels: %d: %s.\n ", returnedCudaError, cudaGetErrorString(returnedCudaError));
			return false;
		}

		return true;
	}

	void reportCudaError(cudaError_t err)
	{
		if (err != cudaSuccess) {
			mexPrintf("CUDA error %d: %s.", err, cudaGetErrorString(err));
            mexErrMsgTxt("ASTRA failed, reboot GPU");
        }
	}



	//
	//float dotproduct3d(cudapitchedptr data, unsigned int x, unsigned int y,
	//                   unsigned int z)
	//{
	//	return astraCUDA3d::dotproduct2d((float*)data.ptr, data.pitch/sizeof(float), x, y*z);
	//}



	int calcNextPowerOfTwo(int _iValue)
	{
		int iOutput = 1;
		while (iOutput < _iValue)
			iOutput *= 2;
		return iOutput;
	}

	double tic()
	{
		return clock();
	}

	double toc(double tstart)
	{
		return (clock() - tstart) / CLOCKS_PER_SEC;
	}



	void printFreeMemory()
	{
		// show memory usage of GPU
		size_t free_byte;
		size_t total_byte;
		cudaError_t cuda_status = cudaMemGetInfo(&free_byte, &total_byte);

		if (cudaSuccess != cuda_status){
			mexPrintf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status));
		}
		double free_db = (double)free_byte;
		double total_db = (double)total_byte;
		double used_db = total_db - free_db;
		mexPrintf("GPU memory usage: used = %g, free = %g MB, total = %g MB\n",
			used_db / 1024.0 / 1024.0, free_db / 1024.0 / 1024.0, total_db / 1024.0 / 1024.0);
	}


	int checkLastError(char * msg)
	{
		cudaError_t cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			char err[512];
			sprintf(err, "astraCUDA3d failed %s: %s. \n", msg, cudaGetErrorString(cudaStatus));
			mexErrMsgTxt(err);
			//mexPrintf(err);
			//mexPrintf("assert \n");
			//ASTRA_CUDA_ASSERT(cudaStatus);
		}
		return 0;
	}

	int dumpArray(char* filename, int width, int height, float *buffer)
	{
		FILE * f;
		int i, j;
		f = fopen(filename, "w");
		for (i = 0; i < height; i++)
		{
			for (j = 0; j < width; j++)
			{
				fprintf(f, "%3.2g\t", buffer[i*width + j]);
				//				fprintf(f, "%i %i\t", i, j);

				//fprintf(f, "%3.2g\t", 1);
			}
			fprintf(f, "\n");
		}
		fclose(f);
		return 0;
	}

	int dumpCudaArray(cudaPitchedPtr Data, int start, int end, char * filename)
	{

		char fname[32], msg[32];
		int width = Data.xsize / sizeof(float);
		int height = Data.ysize;
		int slice_size = width*height*sizeof(float);
		float* buffer = new float[width*height];
		for (int i = start; i < end; i++) {
			cudaMemcpy(buffer, ((float*)Data.ptr) + slice_size*i, slice_size, cudaMemcpyDeviceToHost);
			sprintf(fname, filename, i);
			sprintf(msg, filename, i);
			fprintf(stdout, "%s\n", msg);
			dumpArray(fname, width, height, buffer);
		}
		return 0;
	}

	int writeImageCudaArray(cudaPitchedPtr Data, int start, int end, char * filename)
	{

		char fname[32];
		int width = Data.xsize / sizeof(float);
		int height = Data.ysize;
		int slice_size = width*height*sizeof(float);
		float* buffer = new float[width*height];
		for (int i = start; i < end; i++) {
			cudaMemcpy(buffer, ((float*)Data.ptr) + slice_size*i, slice_size, cudaMemcpyDeviceToHost);
			sprintf(fname, filename, i);
			writeImage(fname, width, height, buffer);
		}
		return 0;
	}

	int writeImage(char * fname, int w, int h, float * data)
	{
		// normalize image 
		float max = 0;
		for (int i = 0; i < w*h; i++)
			if (data[i] > max)
				max = data[i];

		float **x;
		/* allocate the array */
		x = (float **)malloc(h * sizeof *x);
		for (int i = 0; i<h; i++)
			x[i] = (float *)malloc(w * sizeof *x[i]);
		for (int i = 0; i<h; i++)
			for (int j = 0; j < w; j++)
				x[i][j] = data[i*w + j] / max;  // fill the array 

		writeBMPImage(fname, w,h, x,x,x);
		return 0; 
	}



	int writeBMPImage(char * fname, int w, int h, float ** red, float ** green, float ** blue)
	{
	FILE *f;
	unsigned char *img = NULL;
	int filesize = 54 + 3 * w*h;  //w is your image width, h is image height, both int
	if (img)
		free(img);
	img = (unsigned char *)malloc(3 * w*h);
	memset(img, 0, sizeof(img));

	float r, g, b; 
	int x, y; 
	for (int i = 0; i<w; i++)
	{
		for (int j = 0; j<h; j++)
		{
			x = i; y = (h - 1) - j;
			r = red[i][j] * 255;
			g = green[i][j] * 255;
			b = blue[i][j] * 255;
			if (r > 255) r = 255;
			if (g > 255) g = 255;
			if (b > 255) b = 255;
			img[(x + y*w) * 3 + 2] = (unsigned char)(r);
			img[(x + y*w) * 3 + 1] = (unsigned char)(g);
			img[(x + y*w) * 3 + 0] = (unsigned char)(b);
		}
	}

	unsigned char bmpfileheader[14] = { 'B', 'M', 0, 0, 0, 0, 0, 0, 0, 0, 54, 0, 0, 0 };
	unsigned char bmpinfoheader[40] = { 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0 };
	unsigned char bmppad[3] = { 0, 0, 0 };

	bmpfileheader[2] = (unsigned char)(filesize);
	bmpfileheader[3] = (unsigned char)(filesize >> 8);
	bmpfileheader[4] = (unsigned char)(filesize >> 16);
	bmpfileheader[5] = (unsigned char)(filesize >> 24);

	bmpinfoheader[4] = (unsigned char)(w);
	bmpinfoheader[5] = (unsigned char)(w >> 8);
	bmpinfoheader[6] = (unsigned char)(w >> 16);
	bmpinfoheader[7] = (unsigned char)(w >> 24);
	bmpinfoheader[8] = (unsigned char)(h);
	bmpinfoheader[9] = (unsigned char)(h >> 8);
	bmpinfoheader[10] = (unsigned char)(h >> 16);
	bmpinfoheader[11] = (unsigned char)(h >> 24);

	f = fopen(fname, "wb");
	fwrite(bmpfileheader, 1, 14, f);
	fwrite(bmpinfoheader, 1, 40, f);
	for (int i = 0; i < h; i++)
	{
		fwrite(img + (w*(h - i - 1) * 3), 3, w, f);
		fwrite(bmppad, 1, (4 - (w * 3) % 4) % 4, f);
	}
	fclose(f);


	fprintf(stdout, "Saved image %s\n", fname);

	return 0;

	}
}