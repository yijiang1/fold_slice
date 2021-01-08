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
#include <iostream>
#include <list>

#include <cuda.h>
#include "util3d.h"

#include "mex.h"
#include "gpu/mxGPUArray.h"



#ifdef STANDALONE
#include "testutil.h"
#endif

#include "dims3d.h"

typedef texture<float, 3, cudaReadModeElementType> texture3D;

static texture3D gT_par3DVolumeTexture, Xdef0_tex, Ydef0_tex, Zdef0_tex, Xdef1_tex, Ydef1_tex, Zdef1_tex;

#define MAX(x,y) (x>y?x:y);
#define MIN(x,y) (x<y?x:y);


namespace astraCUDA3d {

static const unsigned int g_anglesPerBlock = 4;

// thickness of the slices we're splitting the volume up into
static const unsigned int g_blockSlices = 32;
static const unsigned int g_detBlockU = 32;
static const unsigned int g_detBlockV = 32;

static const unsigned g_MaxAngles = 1024;
__constant__ float gC_RayX[g_MaxAngles];
__constant__ float gC_RayY[g_MaxAngles];
__constant__ float gC_RayZ[g_MaxAngles];
__constant__ float gC_DetSX[g_MaxAngles];
__constant__ float gC_DetSY[g_MaxAngles];
__constant__ float gC_DetSZ[g_MaxAngles];
__constant__ float gC_DetUX[g_MaxAngles];
__constant__ float gC_DetUY[g_MaxAngles];
__constant__ float gC_DetUZ[g_MaxAngles];
__constant__ float gC_DetVX[g_MaxAngles];
__constant__ float gC_DetVY[g_MaxAngles];
__constant__ float gC_DetVZ[g_MaxAngles];
//__constant__ uint8_T gC_use_deform[1];




void __global__ SetVal(float const * const A, float * const B, int const N)
{
	/* Calculate the global linear index, assuming a 1-d grid. */
	int const i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < N) {
		B[i] = A[i];
	}
}




// x=0, y=1, z=2
struct DIR_X {
	__device__ float nSlices(const SDimensions3D& dims) const { return dims.iVolX; }
	__device__ float nDim1(const SDimensions3D& dims) const { return dims.iVolY; }
	__device__ float nDim2(const SDimensions3D& dims) const { return dims.iVolZ; }
	__device__ float c0(float x, float y, float z) const { return x; }
	__device__ float c1(float x, float y, float z) const { return y; }
	__device__ float c2(float x, float y, float z) const { return z; }
	__device__ float tex(float f0, float f1, float f2) const { return tex3D(gT_par3DVolumeTexture, f0, f1, f2); }
	__device__ float texD0x(float f0, float f1, float f2) const { return tex3D(Xdef0_tex, f0, f1, f2); }
	__device__ float texD0y(float f0, float f1, float f2) const { return tex3D(Ydef0_tex, f0, f1, f2); }
	__device__ float texD0z(float f0, float f1, float f2) const { return tex3D(Zdef0_tex, f0, f1, f2); }
	__device__ float texD1x(float f0, float f1, float f2) const { return tex3D(Xdef1_tex, f0, f1, f2); }
	__device__ float texD1y(float f0, float f1, float f2) const { return tex3D(Ydef1_tex, f0, f1, f2); }
	__device__ float texD1z(float f0, float f1, float f2) const { return tex3D(Zdef1_tex, f0, f1, f2); }
	__device__ float x(float f0, float f1, float f2) const { return f0; }
	__device__ float y(float f0, float f1, float f2) const { return f1; }
	__device__ float z(float f0, float f1, float f2) const { return f2; }
	__device__ float offx(const SDimensions3D& dims) const { return dims.iProjU*0.5f; }
	__device__ float offy(const SDimensions3D& dims) const { return dims.iProjV*0.5f; }
};

// y=0, x=1, z=2
struct DIR_Y {
	__device__ float nSlices(const SDimensions3D& dims) const { return dims.iVolY; }
	__device__ float nDim1(const SDimensions3D& dims) const { return dims.iVolX; }
	__device__ float nDim2(const SDimensions3D& dims) const { return dims.iVolZ; }
	__device__ float c0(float x, float y, float z) const { return y; }
	__device__ float c1(float x, float y, float z) const { return x; }
	__device__ float c2(float x, float y, float z) const { return z; }
	__device__ float tex(float f0, float f1, float f2) const { return tex3D(gT_par3DVolumeTexture, f1, f0, f2); }
	__device__ float texD0x(float f0, float f1, float f2) const { return tex3D(Ydef0_tex, f1, f0, f2); }
	__device__ float texD0y(float f0, float f1, float f2) const { return tex3D(Xdef0_tex, f1, f0, f2); }
	__device__ float texD0z(float f0, float f1, float f2) const { return tex3D(Zdef0_tex, f1, f0, f2); }
	__device__ float texD1x(float f0, float f1, float f2) const { return tex3D(Ydef1_tex, f1, f0, f2); }
	__device__ float texD1y(float f0, float f1, float f2) const { return tex3D(Xdef1_tex, f1, f0, f2); }
	__device__ float texD1z(float f0, float f1, float f2) const { return tex3D(Zdef1_tex, f1, f0, f2); }
	__device__ float x(float f0, float f1, float f2) const { return f1; }
	__device__ float y(float f0, float f1, float f2) const { return f0; }
	__device__ float z(float f0, float f1, float f2) const { return f2; }
	__device__ float offx(const SDimensions3D& dims) const { return dims.iProjU*0.5f; }
	__device__ float offy(const SDimensions3D& dims) const { return dims.iProjV*0.5f; }
};

// z=0, x=1, y=2
struct DIR_Z {
	__device__ float nSlices(const SDimensions3D& dims) const { return dims.iVolZ; }
	__device__ float nDim1(const SDimensions3D& dims) const { return dims.iVolX; }
	__device__ float nDim2(const SDimensions3D& dims) const { return dims.iVolY; }
	__device__ float c0(float x, float y, float z) const { return z; }
	__device__ float c1(float x, float y, float z) const { return x; }
	__device__ float c2(float x, float y, float z) const { return y; }
	__device__ float tex(float f0, float f1, float f2) const { return tex3D(gT_par3DVolumeTexture, f1, f2, f0); }
	__device__ float texD0x(float f0, float f1, float f2) const { return tex3D(Zdef0_tex, f1, f2, f0); }
	__device__ float texD0y(float f0, float f1, float f2) const { return tex3D(Xdef0_tex, f1, f2, f0); }
	__device__ float texD0z(float f0, float f1, float f2) const { return tex3D(Ydef0_tex, f1, f2, f0); }
	__device__ float texD1x(float f0, float f1, float f2) const { return tex3D(Zdef1_tex, f1, f2, f0); }
	__device__ float texD1y(float f0, float f1, float f2) const { return tex3D(Xdef1_tex, f1, f2, f0); }
	__device__ float texD1z(float f0, float f1, float f2) const { return tex3D(Ydef1_tex, f1, f2, f0); }
	__device__ float x(float f0, float f1, float f2) const { return f1; }
	__device__ float y(float f0, float f1, float f2) const { return f2; }
	__device__ float z(float f0, float f1, float f2) const { return f0; }
	__device__ float offx(const SDimensions3D& dims) const { return dims.iProjU*0.5f; }
	__device__ float offy(const SDimensions3D& dims) const { return dims.iProjV*0.5f; }
};



// threadIdx: x = u detector
//            y = relative angle
// blockIdx:  x = u/v detector
//            y = angle block





template<class COORD>
__global__ void par3D_FP_t(float* D_projData, unsigned int projPitch,
                           unsigned int startSlice,
                           unsigned int startAngle, unsigned int endAngle,
                           const SDimensions3D dims, float fOutputScale, const bool use_deform, const bool linear_deform_model)
{
	COORD c;


	int angle = startAngle + blockIdx.y * g_anglesPerBlock + threadIdx.y;
	if (angle >= endAngle)
		return;



	const float fRayX = gC_RayX[angle];
	const float fRayY = gC_RayY[angle];
	const float fRayZ = gC_RayZ[angle];
	const float fDetUX = gC_DetUX[angle];
	const float fDetUY = gC_DetUY[angle];
	const float fDetUZ = gC_DetUZ[angle];
	const float fDetVX = gC_DetVX[angle];
	const float fDetVY = gC_DetVY[angle];
	const float fDetVZ = gC_DetVZ[angle];
	const float fDetSX = gC_DetSX[angle] + 0.5f * fDetUX + 0.5f * fDetVX;
	const float fDetSY = gC_DetSY[angle] + 0.5f * fDetUY + 0.5f * fDetVY;
	const float fDetSZ = gC_DetSZ[angle] + 0.5f * fDetUZ + 0.5f * fDetVZ;



	if (c.c0(fRayX, fRayY, fRayZ) == 0)
		return;

	const int detectorU = (blockIdx.x%((dims.iProjU+g_detBlockU-1)/g_detBlockU)) * g_detBlockU + threadIdx.x;

    if (detectorU >= dims.iProjU) 
        return; 
    

	const int startDetectorV = (blockIdx.x/((dims.iProjU+g_detBlockU-1)/g_detBlockU)) * g_detBlockV;
	int endDetectorV = startDetectorV + g_detBlockV;
	if (endDetectorV > dims.iProjV)
		endDetectorV = dims.iProjV;

	int endSlice = startSlice + g_blockSlices;
	if (endSlice > c.nSlices(dims))
		endSlice = c.nSlices(dims);

	// FIXME 
	/*if (endSlice < startSlice - 1)
		return;*/

    float angle_ratio = (float)angle / (float)dims.iProjAngles ; 


	for (int detectorV = startDetectorV; detectorV < endDetectorV; ++detectorV)
	{
		/* Trace ray in direction Ray to (detectorU,detectorV) from  */
		/* X = startSlice to X = endSlice                            */

		const float fDetX = fDetSX + (detectorU*fDetUX + detectorV*fDetVX);
		const float fDetY = fDetSY + (detectorU*fDetUY + detectorV*fDetVY);
		const float fDetZ = fDetSZ + (detectorU*fDetUZ + detectorV*fDetVZ);

		/*        (x)   ( 1)       ( 0)    */
		/* ray:   (y) = (ay) * x + (by)    */
		/*        (z)   (az)       (bz)    */

		const float a1 = c.c1(fRayX,fRayY,fRayZ) / c.c0(fRayX,fRayY,fRayZ);
		const float a2 = c.c2(fRayX,fRayY,fRayZ) / c.c0(fRayX,fRayY,fRayZ);
		const float b1 = c.c1(fDetX,fDetY,fDetZ) - a1 * c.c0(fDetX,fDetY,fDetZ);
		const float b2 = c.c2(fDetX,fDetY,fDetZ) - a2 * c.c0(fDetX,fDetY,fDetZ);

		const float fDistCorr = sqrt(a1*a1+a2*a2+1.0f) * fOutputScale;

		float fVal = 0.0f;

		//float f0 =  startSlice + 0.5f;
		//float f1 =  a1 * (startSlice - 0.5f*c.nSlices(dims) + 0.5f) + b1 + 0.5f*c.nDim1(dims) - 0.5f + 0.5f;
		//float f2 =  a2 * (startSlice - 0.5f*c.nSlices(dims) + 0.5f) + b2 + 0.5f*c.nDim2(dims) - 0.5f + 0.5f;

        bool is_inside; 
        int lim0, lim1, lim2; 
        lim0 = c.nSlices(dims); 
        lim1 = c.nDim1(dims); 
        lim2 = c.nDim2(dims);
        const float offset = 0.5*lim0; 

        // calculate minimal distance needed to get the subprojection, important for laminography and large projection size 

		int startSlice_tmp = startSlice; 
		int endSlice_tmp =  endSlice;

		if (a1 > 0)
		{
            startSlice_tmp = MAX(startSlice_tmp, floor((-0.5*lim1-b1-0.5f)/a1+offset-1.0f)); 
              endSlice_tmp = MIN(endSlice_tmp,    ceil((+0.5*lim1-b1+0.5f)/a1+offset+1.0f)); 
		}
		else if (a1 < 0)
		{
            startSlice_tmp = MAX(startSlice_tmp, floor((+0.5*lim1-b1+0.5f)/a1+offset-1.0f)); 
              endSlice_tmp = MIN(endSlice_tmp,    ceil((-0.5*lim1-b1-0.5f)/a1+offset+1.0f)); 
		}
		if (a2 > 0)
		{
           startSlice_tmp = MAX(startSlice_tmp, floor((-0.5*lim2-b2-0.5f)/a2+offset-1.0f));
             endSlice_tmp = MIN(endSlice_tmp,    ceil((+0.5*lim2-b2+0.5f)/a2+offset+1.0f));  
  		}
	    else if (a2 < 0)
		{
          startSlice_tmp = MAX(startSlice_tmp,  floor((+0.5*lim2-b2+0.5f)/a2+offset-1.0f));
            endSlice_tmp = MIN(endSlice_tmp,     ceil((-0.5*lim2-b2-0.5f)/a2+offset+1.0f));  
		}

        endSlice_tmp   = MIN(endSlice_tmp,     endSlice); 
        endSlice_tmp   = MAX(endSlice_tmp,     0); 

        startSlice_tmp = MAX(startSlice_tmp, startSlice); 
        startSlice_tmp = MIN(startSlice_tmp, endSlice_tmp);



		float f0 =  startSlice_tmp + 0.5f;
		float f1 =  a1 * (startSlice_tmp - offset+0.5f) + b1 + 0.5f*c.nDim1(dims);
		float f2 =  a2 * (startSlice_tmp - offset+0.5f) + b2 + 0.5f*c.nDim2(dims);


        float f0s, f1s, f2s;  // shifted coordinates 
        float f0n, f1n, f2n; // normalized coordinates 

		// 87% of the execution time 
		for (int s = startSlice_tmp; s < endSlice_tmp; ++s)
		{
        if (use_deform)  {
            f0n = f0/lim0; // normalized coordinates 
            f1n = f1/lim1; 
            f2n = f2/lim2; 
            
            // load deformed coordinates 
            if (!linear_deform_model) {
                f0s =  f0 - c.texD0x(f0n, f1n, f2n); 
                f1s =  f1 - c.texD0y(f0n, f1n, f2n);
                f2s =  f2 - c.texD0z(f0n, f1n, f2n);
            } else {
                f0s =  f0 - (c.texD0x(f0n, f1n, f2n) * (1-angle_ratio) + (angle_ratio)*c.texD1x(f0n, f1n, f2n)); 
                f1s =  f1 - (c.texD0y(f0n, f1n, f2n) * (1-angle_ratio) + (angle_ratio)*c.texD1y(f0n, f1n, f2n)); 
                f2s =  f2 - (c.texD0z(f0n, f1n, f2n) * (1-angle_ratio) + (angle_ratio)*c.texD1z(f0n, f1n, f2n)); 
            }
            // get trilinear interpolation in the shifted coordinates 
            fVal += c.tex(f0s, f1s, f2s); 
        } else {

             is_inside = (f0 > 0 && f1 > 0 && f2 > 0 && f0 < lim0 && f1 < lim1 && f2 < lim2 ); 
            //  fVal += (is_inside ? c.tex(f0, f1, f2) : 0); // skip textures on boundaries 
            //fVal += c.tex(f0, f1, f2) == 0; 
            //fVal += is_inside == 0; 

            fVal += c.tex(f0, f1, f2); // fastest seems to be let texture memory to handle boundaries 
        }

            // move to the next pixel 
			f0 += 1.0f;
			f1 += a1;
			f2 += a2;
		}

		fVal *= fDistCorr;

		// !! 10% of the execution time 
		//D_projData[(detectorV*dims.iProjAngles + angle)*projPitch + detectorU] += fVal;
		atomicAdd(&D_projData[(detectorV*dims.iProjAngles + angle)*projPitch + detectorU], fVal);

	}
}

// Supersampling version
template<class COORD>
__global__ void par3D_FP_SS_t(float* D_projData, unsigned int projPitch,
                              unsigned int startSlice,
                              unsigned int startAngle, unsigned int endAngle,
                              const SDimensions3D dims, float fOutputScale)
{
	COORD c;

	int angle = startAngle + blockIdx.y * g_anglesPerBlock + threadIdx.y;
	if (angle >= endAngle)
		return;

	const float fRayX = gC_RayX[angle];
	const float fRayY = gC_RayY[angle];
	const float fRayZ = gC_RayZ[angle];
	const float fDetUX = gC_DetUX[angle];
	const float fDetUY = gC_DetUY[angle];
	const float fDetUZ = gC_DetUZ[angle];
	const float fDetVX = gC_DetVX[angle];
	const float fDetVY = gC_DetVY[angle];
	const float fDetVZ = gC_DetVZ[angle];
	const float fDetSX = gC_DetSX[angle] + 0.5f * fDetUX + 0.5f * fDetVX;
	const float fDetSY = gC_DetSY[angle] + 0.5f * fDetUY + 0.5f * fDetVY;
	const float fDetSZ = gC_DetSZ[angle] + 0.5f * fDetUZ + 0.5f * fDetVZ;



	const int detectorU = (blockIdx.x%((dims.iProjU+g_detBlockU-1)/g_detBlockU)) * g_detBlockU + threadIdx.x;
	const int startDetectorV = (blockIdx.x/((dims.iProjU+g_detBlockU-1)/g_detBlockU)) * g_detBlockV;
	int endDetectorV = startDetectorV + g_detBlockV;
	if (endDetectorV > dims.iProjV)
		endDetectorV = dims.iProjV;

	int endSlice = startSlice + g_blockSlices;
	if (endSlice > c.nSlices(dims))
		endSlice = c.nSlices(dims);

	const float fSubStep = 1.0f/dims.iRaysPerDetDim;

	for (int detectorV = startDetectorV; detectorV < endDetectorV; ++detectorV)
	{

		float fV = 0.0f;

		float fdU = detectorU - 0.5f + 0.5f*fSubStep;
		for (int iSubU = 0; iSubU < dims.iRaysPerDetDim; ++iSubU, fdU+=fSubStep) {
		float fdV = detectorV - 0.5f + 0.5f*fSubStep;
		for (int iSubV = 0; iSubV < dims.iRaysPerDetDim; ++iSubV, fdV+=fSubStep) {

		/* Trace ray in direction Ray to (detectorU,detectorV) from  */
		/* X = startSlice to X = endSlice                            */

		const float fDetX = fDetSX + fdU*fDetUX + fdV*fDetVX;
		const float fDetY = fDetSY + fdU*fDetUY + fdV*fDetVY;
		const float fDetZ = fDetSZ + fdU*fDetUZ + fdV*fDetVZ;

		/*        (x)   ( 1)       ( 0)    */
		/* ray:   (y) = (ay) * x + (by)    */
		/*        (z)   (az)       (bz)    */


		const float a1 = c.c1(fRayX,fRayY,fRayZ) / c.c0(fRayX,fRayY,fRayZ);
		const float a2 = c.c2(fRayX,fRayY,fRayZ) / c.c0(fRayX,fRayY,fRayZ);
		const float b1 = c.c1(fDetX,fDetY,fDetZ) - a1 * c.c0(fDetX,fDetY,fDetZ);
		const float b2 = c.c2(fDetX,fDetY,fDetZ) - a2 * c.c0(fDetX,fDetY,fDetZ);

		const float fDistCorr = sqrt(a1*a1+a2*a2+1.0f) * fOutputScale;

		float fVal = 0.0f;

		float f0 = startSlice + 0.5f;
		float f1 = a1 * (startSlice - 0.5f*c.nSlices(dims) + 0.5f) + b1 + 0.5f*c.nDim1(dims) - 0.5f + 0.5f;
		float f2 = a2 * (startSlice - 0.5f*c.nSlices(dims) + 0.5f) + b2 + 0.5f*c.nDim2(dims) - 0.5f + 0.5f;



		for (int s = startSlice; s < endSlice; ++s)
		{
			fVal += c.tex(f0, f1, f2);
			f0 += 1.0f;
			f1 += a1 ;
			// f2 += a2;
		}

		fVal *= fDistCorr;
		fV += fVal;

		}
		}

		D_projData[(detectorV*dims.iProjAngles+angle)*projPitch+detectorU] += fV / (dims.iRaysPerDetDim * dims.iRaysPerDetDim);
	}
}


__device__ float dirWeights(float fX, float fN) {
	if (fX <= -0.5f) // outside image on left
		return 0.0f;
	if (fX <= 0.5f) // half outside image on left
		return (fX + 0.5f) * (fX + 0.5f);
	if (fX <= fN - 0.5f) { // inside image
		float t = fX + 0.5f - floorf(fX + 0.5f);
		return 1; //   t*t + (1 - t)*(1 - t);
	}
	if (fX <= fN + 0.5f) // half outside image on right
		return (fN + 0.5f - fX) * (fN + 0.5f - fX);
	return 0.0f; // outside image on right
}

template<class COORD>
__global__ void par3D_FP_SumSqW_t(float* D_projData, unsigned int projPitch,
                                  unsigned int startSlice,
                                  unsigned int startAngle, unsigned int endAngle,
                                  const SDimensions3D dims, float fOutputScale)
{
	COORD c;

	int angle = startAngle + blockIdx.y * g_anglesPerBlock + threadIdx.y;
	if (angle >= endAngle)
		return;

	const float fRayX = gC_RayX[angle];
	const float fRayY = gC_RayY[angle];
	const float fRayZ = gC_RayZ[angle];
	const float fDetUX = gC_DetUX[angle];
	const float fDetUY = gC_DetUY[angle];
	const float fDetUZ = gC_DetUZ[angle];
	const float fDetVX = gC_DetVX[angle];
	const float fDetVY = gC_DetVY[angle];
	const float fDetVZ = gC_DetVZ[angle];
	const float fDetSX = gC_DetSX[angle] + 0.5f * fDetUX + 0.5f * fDetVX;
	const float fDetSY = gC_DetSY[angle] + 0.5f * fDetUY + 0.5f * fDetVY;
	const float fDetSZ = gC_DetSZ[angle] + 0.5f * fDetUZ + 0.5f * fDetVZ;



	const int detectorU = (blockIdx.x%((dims.iProjU+g_detBlockU-1)/g_detBlockU)) * g_detBlockU + threadIdx.x;
	const int startDetectorV = (blockIdx.x/((dims.iProjU+g_detBlockU-1)/g_detBlockU)) * g_detBlockV;
	int endDetectorV = startDetectorV + g_detBlockV;
	if (endDetectorV > dims.iProjV)
		endDetectorV = dims.iProjV;

	int endSlice = startSlice + g_blockSlices;
	if (endSlice > c.nSlices(dims))
		endSlice = c.nSlices(dims);

	for (int detectorV = startDetectorV; detectorV < endDetectorV; ++detectorV)
	{
		/* Trace ray in direction Ray to (detectorU,detectorV) from  */
		/* X = startSlice to X = endSlice                            */

		const float fDetX = fDetSX + detectorU*fDetUX + detectorV*fDetVX;
		const float fDetY = fDetSY + detectorU*fDetUY + detectorV*fDetVY;
		const float fDetZ = fDetSZ + detectorU*fDetUZ + detectorV*fDetVZ;

		/*        (x)   ( 1)       ( 0)    */
		/* ray:   (y) = (ay) * x + (by)    */
		/*        (z)   (az)       (bz)    */

		const float a1 = c.c1(fRayX,fRayY,fRayZ) / c.c0(fRayX,fRayY,fRayZ);
		const float a2 = c.c2(fRayX,fRayY,fRayZ) / c.c0(fRayX,fRayY,fRayZ);
		const float b1 = c.c1(fDetX,fDetY,fDetZ) - a1 * c.c0(fDetX,fDetY,fDetZ);
		const float b2 = c.c2(fDetX,fDetY,fDetZ) - a2 * c.c0(fDetX,fDetY,fDetZ);

		const float fDistCorr = sqrt(a1*a1+a2*a2+1.0f) * fOutputScale;

		float fVal = 0.0f;

		float f0 = startSlice + 0.5f;
		float f1 = a1 * (startSlice - 0.5f*c.nSlices(dims) + 0.5f) + b1 + 0.5f*c.nDim1(dims) - 0.5f + 0.5f;
		float f2 = a2 * (startSlice - 0.5f*c.nSlices(dims) + 0.5f) + b2 + 0.5f*c.nDim2(dims) - 0.5f + 0.5f;

		for (int s = startSlice; s < endSlice; ++s)
		{
			fVal += dirWeights(f1, c.nDim1(dims)) * dirWeights(f2, c.nDim2(dims)) * fDistCorr * fDistCorr;
			f0 += 1.0f;
			f1 += a1;
			f2 += a2;
		}

		D_projData[(detectorV*dims.iProjAngles+angle)*projPitch+detectorU] += fVal;
	}
}

// Supersampling version
// TODO


bool Par3DFP_Array_internal(cudaPitchedPtr D_projData,
                   const SDimensions3D& dims, unsigned int angleCount, const SPar3DProjection* angles,
                   float fOutputScale, const bool use_deform, const bool linear_deform_model)
{



	// transfer angles to constant memory
	float* tmp = new float[dims.iProjAngles];

#define TRANSFER_TO_CONSTANT(name) do { for (unsigned int i = 0; i < angleCount; ++i) tmp[i] = (float)angles[i].f##name ; cudaMemcpyToSymbol(gC_##name, tmp, angleCount*sizeof(float), 0, cudaMemcpyHostToDevice); } while (0)

	TRANSFER_TO_CONSTANT(RayX);
	TRANSFER_TO_CONSTANT(RayY);
	TRANSFER_TO_CONSTANT(RayZ);
	TRANSFER_TO_CONSTANT(DetSX);
	TRANSFER_TO_CONSTANT(DetSY);
	TRANSFER_TO_CONSTANT(DetSZ);
	TRANSFER_TO_CONSTANT(DetUX);
	TRANSFER_TO_CONSTANT(DetUY);
	TRANSFER_TO_CONSTANT(DetUZ);
	TRANSFER_TO_CONSTANT(DetVX);
	TRANSFER_TO_CONSTANT(DetVY);
	TRANSFER_TO_CONSTANT(DetVZ);

#undef TRANSFER_TO_CONSTANT

	delete[] tmp;

	std::list<cudaStream_t> streams;
	dim3 dimBlock(g_detBlockU, g_anglesPerBlock); // region size, angles

	// Run over all angles, grouping them into groups of the same
	// orientation (roughly horizontal vs. roughly vertical).
	// Start a stream of grids for each such group.

	unsigned int blockStart = 0;
	unsigned int blockEnd = 0;
	int blockDirection = 0;


	for (unsigned int a = 0; a <= angleCount; ++a) {
		int dir = -1;
		if (a != dims.iProjAngles) {
			float dX = fabsf(angles[a].fRayX);
			float dY = fabsf(angles[a].fRayY);
			float dZ = fabsf(angles[a].fRayZ);

			if (dX >= dY && dX >= dZ)
				dir = 0;
			else if (dY >= dX && dY >= dZ)
				dir = 1;
			else
				dir = 2;
		}

		if (a == angleCount || dir != blockDirection) {
			// block done

			blockEnd = a;
			if (blockStart != blockEnd) {

				dim3 dimGrid(
				             ((dims.iProjU+g_detBlockU-1)/g_detBlockU)*((dims.iProjV+g_detBlockV-1)/g_detBlockV),
                                    (blockEnd-blockStart+g_anglesPerBlock-1)/g_anglesPerBlock);
				// TODO: check if we can't immediately
				//       destroy the stream after use


				cudaStream_t stream;
				cudaStreamCreate(&stream);
				streams.push_back(stream);

			    //mexPrintf("angle block: %d to %d, %d (%dx%d, %dx%d)\n", blockStart, blockEnd, blockDirection, dimGrid.x, dimGrid.y, dimBlock.x, dimBlock.y);
				//mexPrintf(" Nelements %i ", (dims.iProjU)*(dims.iProjV)*(dims.iProjAngles)); 
			

				if (blockDirection == 0) {
					for (unsigned int i = 0; i < dims.iVolX; i += g_blockSlices)
						if (dims.iRaysPerDetDim == 1)
							par3D_FP_t<DIR_X><<<dimGrid, dimBlock, 0, stream>>>((float*)D_projData.ptr, D_projData.pitch/sizeof(float), i, blockStart, blockEnd, dims, fOutputScale, use_deform, linear_deform_model);
						else
							par3D_FP_SS_t<DIR_X><<<dimGrid, dimBlock, 0, stream>>>((float*)D_projData.ptr, D_projData.pitch/sizeof(float), i, blockStart, blockEnd, dims, fOutputScale);
				} else if (blockDirection == 1) {
					for (unsigned int i = 0; i < dims.iVolY; i += g_blockSlices)
						if (dims.iRaysPerDetDim == 1)
                            par3D_FP_t<DIR_Y><<<dimGrid, dimBlock, 0, stream>>>((float*)D_projData.ptr, D_projData.pitch/sizeof(float), i, blockStart, blockEnd, dims, fOutputScale, use_deform, linear_deform_model);
						else
							par3D_FP_SS_t<DIR_Y><<<dimGrid, dimBlock, 0, stream>>>((float*)D_projData.ptr, D_projData.pitch/sizeof(float), i, blockStart, blockEnd, dims, fOutputScale);
				} else if (blockDirection == 2) {
					for (unsigned int i = 0; i < dims.iVolZ; i += g_blockSlices)
                        if (dims.iRaysPerDetDim == 1)
						par3D_FP_t<DIR_Z><<<dimGrid, dimBlock, 0, stream>>>((float*)D_projData.ptr, D_projData.pitch/sizeof(float), i, blockStart, blockEnd, dims, fOutputScale, use_deform, linear_deform_model);
						else
							par3D_FP_SS_t<DIR_Z><<<dimGrid, dimBlock, 0, stream>>>((float*)D_projData.ptr, D_projData.pitch/sizeof(float), i, blockStart, blockEnd, dims, fOutputScale);
				}

			}

			blockDirection = dir;
			blockStart = a;
		}
	}

	cudaThreadSynchronize(); 

	for (std::list<cudaStream_t>::iterator iter = streams.begin(); iter != streams.end(); ++iter)
		cudaStreamDestroy(*iter);
	

	streams.clear();

	cudaTextForceKernelsCompletion();

	return true;
}

bool Par3DFP(cudaPitchedPtr D_volumeData,
             cudaPitchedPtr D_projData,
             const SDimensions3D& dims, const SPar3DProjection* angles,
             float fOutputScale, DeformField DF)
{



	checkLastError("before allocateVolumeArray");
	/*printFreeMemory(); 
	mexPrintf("Allocate memory\n");*/
	// transfer volume to array


	if (dims.iVolX*dims.iVolY*dims.iVolZ * 4 > 1024e6)
	{
		mexPrintf("Volume exceeded maximal size of texture 1024MB \n");
		return 1;
	}

	cudaArray* cuArray = allocateVolumeArray(dims);
	//mexPrintf("Allocate memory done\n"); 
	//printFreeMemory();

	checkLastError("after allocateVolumeArray");

	//mexPrintf("transferVolumeToArray\n"); 

    transferVolumeToArray(D_volumeData, cuArray, dims);


	checkLastError("after transferVolumeToArray\n \n ");
	//printFreeMemory();

	bindDataTexture(cuArray, gT_par3DVolumeTexture,cudaAddressModeBorder, false);

	//mexPrintf("bindDataTexture done \n");


	checkLastError("after bindDataTexture");
	//printFreeMemory();

	//mexPrintf("preoparation finieshe \n");
    cudaArray * cuArrX0, *cuArrY0, *cuArrZ0, *cuArrX1, *cuArrY1, *cuArrZ1 ; 


    if (DF.use_deform) {
       // mexPrintf("transferDeformationToArray\n");

        cuArrX0 = transferDeformationToArray(DF.X0);
        cuArrY0 = transferDeformationToArray(DF.Y0);
        cuArrZ0 = transferDeformationToArray(DF.Z0);
        bindDataTexture(cuArrX0, Xdef0_tex,cudaAddressModeClamp, true);
        bindDataTexture(cuArrY0, Ydef0_tex,cudaAddressModeClamp, true);
        bindDataTexture(cuArrZ0, Zdef0_tex,cudaAddressModeClamp, true);
        if (DF.use_linear_model) {
            cuArrX1 = transferDeformationToArray(DF.X1);
            cuArrY1 = transferDeformationToArray(DF.Y1);
            cuArrZ1 = transferDeformationToArray(DF.Z1);
            bindDataTexture(cuArrX1, Xdef1_tex,cudaAddressModeClamp, true);
            bindDataTexture(cuArrY1, Ydef1_tex,cudaAddressModeClamp, true);
            bindDataTexture(cuArrZ1, Zdef1_tex,cudaAddressModeClamp, true);
        }
    }


	bool ret;

    // ONLY A LIMITED RANGE OF ANGLES IS AVAILIBLE INSIDE !!!!!
    checkLastError("before allocateVolumeArray");


	// 97% of time spent in Par3DFP_Array_internal 
		ret = Par3DFP_Array_internal(D_projData,
		                             dims, dims.iProjAngles, angles,
		                             fOutputScale, DF.use_deform, DF.use_linear_model);
    checkLastError("after allocateVolumeArray");

	cudaFreeArray(cuArray);
	checkLastError("after cudaFreeArray");

	// THIS WAS BUG IN ASTRA !!!! 
	cudaUnbindTexture(gT_par3DVolumeTexture);
	checkLastError("cudaUnbindTexture");

    if (DF.use_deform) {
        cudaFreeArray(cuArrX0);
        cudaFreeArray(cuArrY0);
        cudaFreeArray(cuArrZ0);
        cudaUnbindTexture(Xdef0_tex);
        cudaUnbindTexture(Ydef0_tex);
        cudaUnbindTexture(Zdef0_tex);
        if (DF.use_linear_model) {
            cudaFreeArray(cuArrX1);
            cudaFreeArray(cuArrY1);
            cudaFreeArray(cuArrZ1);
            cudaUnbindTexture(Xdef1_tex);
            cudaUnbindTexture(Ydef1_tex);
            cudaUnbindTexture(Zdef1_tex);
        }
        checkLastError("unbind deforms");
    }






	return ret;
}



bool Par3DFP_SumSqW(cudaPitchedPtr D_volumeData,
                    cudaPitchedPtr D_projData,
                    const SDimensions3D& dims, const SPar3DProjection* angles,
                    float fOutputScale)
{
	// transfer angles to constant memory
	float* tmp = new float[dims.iProjAngles];

#define TRANSFER_TO_CONSTANT(name) do { for (unsigned int i = 0; i < dims.iProjAngles; ++i) tmp[i] = angles[i].f##name ; cudaMemcpyToSymbol(gC_##name, tmp, dims.iProjAngles*sizeof(float), 0, cudaMemcpyHostToDevice); } while (0)

	TRANSFER_TO_CONSTANT(RayX);
	TRANSFER_TO_CONSTANT(RayY);
	TRANSFER_TO_CONSTANT(RayZ);
	TRANSFER_TO_CONSTANT(DetSX);
	TRANSFER_TO_CONSTANT(DetSY);
	TRANSFER_TO_CONSTANT(DetSZ);
	TRANSFER_TO_CONSTANT(DetUX);
	TRANSFER_TO_CONSTANT(DetUY);
	TRANSFER_TO_CONSTANT(DetUZ);
	TRANSFER_TO_CONSTANT(DetVX);
	TRANSFER_TO_CONSTANT(DetVY);
	TRANSFER_TO_CONSTANT(DetVZ);

#undef TRANSFER_TO_CONSTANT

	delete[] tmp;

	std::list<cudaStream_t> streams;
	dim3 dimBlock(g_detBlockU, g_anglesPerBlock); // region size, angles

	// Run over all angles, grouping them into groups of the same
	// orientation (roughly horizontal vs. roughly vertical).
	// Start a stream of grids for each such group.

	unsigned int blockStart = 0;
	unsigned int blockEnd = 0;
	int blockDirection = 0;

	// timeval t;
	// tic(t);

	for (unsigned int a = 0; a <= dims.iProjAngles; ++a) {
		int dir;
		if (a != dims.iProjAngles) {
			float dX = fabsf(angles[a].fRayX);
			float dY = fabsf(angles[a].fRayY);
			float dZ = fabsf(angles[a].fRayZ);

			if (dX >= dY && dX >= dZ)
				dir = 0;
			else if (dY >= dX && dY >= dZ)
				dir = 1;
			else
				dir = 2;
		}

		if (a == dims.iProjAngles || dir != blockDirection) {
			// block done

			blockEnd = a;
			if (blockStart != blockEnd) {

				dim3 dimGrid(
				             ((dims.iProjU+g_detBlockU-1)/g_detBlockU)*((dims.iProjV+g_detBlockV-1)/g_detBlockV),
(blockEnd-blockStart+g_anglesPerBlock-1)/g_anglesPerBlock);
				// TODO: check if we can't immediately
				//       destroy the stream after use
				cudaStream_t stream;
				cudaStreamCreate(&stream);
				streams.push_back(stream);

				//printf("angle block: %d to %d, %d (%dx%d, %dx%d)\n", blockStart, blockEnd, blockDirection, dimGrid.x, dimGrid.y, dimBlock.x, dimBlock.y);

				if (blockDirection == 0) {
					for (unsigned int i = 0; i < dims.iVolX; i += g_blockSlices)
						if (dims.iRaysPerDetDim == 1)
							par3D_FP_SumSqW_t<DIR_X><<<dimGrid, dimBlock, 0, stream>>>((float*)D_projData.ptr, D_projData.pitch/sizeof(float), i, blockStart, blockEnd, dims, fOutputScale);
						else
#if 0
							par3D_FP_SS_SumSqW_dirX<<<dimGrid, dimBlock, 0, stream>>>((float*)D_projData.ptr, D_projData.pitch/sizeof(float), i, blockStart, blockEnd, dims, fOutputScale);
#else
							assert(false);
#endif
				} else if (blockDirection == 1) {
					for (unsigned int i = 0; i < dims.iVolY; i += g_blockSlices)
						if (dims.iRaysPerDetDim == 1)
							par3D_FP_SumSqW_t<DIR_Y><<<dimGrid, dimBlock, 0, stream>>>((float*)D_projData.ptr, D_projData.pitch/sizeof(float), i, blockStart, blockEnd, dims, fOutputScale);
						else
#if 0
							par3D_FP_SS_SumSqW_dirY<<<dimGrid, dimBlock, 0, stream>>>((float*)D_projData.ptr, D_projData.pitch/sizeof(float), i, blockStart, blockEnd, dims, fOutputScale);
#else
							assert(false);
#endif
				} else if (blockDirection == 2) {
					for (unsigned int i = 0; i < dims.iVolZ; i += g_blockSlices)
						if (dims.iRaysPerDetDim == 1)
							par3D_FP_SumSqW_t<DIR_Z><<<dimGrid, dimBlock, 0, stream>>>((float*)D_projData.ptr, D_projData.pitch/sizeof(float), i, blockStart, blockEnd, dims, fOutputScale);
						else
#if 0
							par3D_FP_SS_SumSqW_dirZ<<<dimGrid, dimBlock, 0, stream>>>((float*)D_projData.ptr, D_projData.pitch/sizeof(float), i, blockStart, blockEnd, dims, fOutputScale);
#else
							assert(false);
#endif
				}

			}

			blockDirection = dir;
			blockStart = a;
		}
	}

	for (std::list<cudaStream_t>::iterator iter = streams.begin(); iter != streams.end(); ++iter)
		cudaStreamDestroy(*iter);

	streams.clear();

	cudaTextForceKernelsCompletion();


	// printf("%f\n", toc(t));

	return true;
}

}

