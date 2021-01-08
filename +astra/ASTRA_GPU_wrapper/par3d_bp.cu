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

#ifdef STANDALONE
#include "par3d_fp.h"
#include "testutil.h"
#endif

#include "dims3d.h"

typedef texture<float, 3, cudaReadModeElementType> texture3D;

static texture3D gT_par3DProjTexture, Xdef0_tex, Ydef0_tex, Zdef0_tex, Xdef1_tex, Ydef1_tex, Zdef1_tex;

namespace astraCUDA3d {

#define ZSIZE 6
static const unsigned int g_volBlockZ = ZSIZE;

static const unsigned int g_anglesPerBlock = 32;
static const unsigned int g_volBlockX = 16;
static const unsigned int g_volBlockY = 32;

static const unsigned g_MaxAngles = 1024;

__constant__ float gC_C[8*g_MaxAngles];

#define MAX(x,y) (x>y?x:y);
#define MIN(x,y) (x<y?x:y);
#define ABS(x) (x>0?x:-x);






__global__ void dev_par3D_BP(void* D_volData, unsigned int volPitch, 
        int startAngle, int angleOffset, const SDimensions3D dims, 
        float fOutputScale, bool use_deform, bool linear_deform_model)
{
	float* volData = (float*)D_volData;

	int endAngle = startAngle + g_anglesPerBlock;
	if (endAngle > dims.iProjAngles - angleOffset)
		endAngle = dims.iProjAngles - angleOffset;

	// threadIdx: x = rel x
	//            y = rel y

	// blockIdx:  x = x + y
	//            y = z


	const int X = blockIdx.x % ((dims.iVolX+g_volBlockX-1)/g_volBlockX) * g_volBlockX + threadIdx.x;
	const int Y = blockIdx.x / ((dims.iVolX+g_volBlockX-1)/g_volBlockX) * g_volBlockY + threadIdx.y;

	if (X >= dims.iVolX)
		return;
	if (Y >= dims.iVolY)
		return;

	const int startZ = blockIdx.y * g_volBlockZ;

    const float limX = dims.iVolX; 
    const float limY = dims.iVolY; 
    const float limZ = dims.iVolZ; 

	float fX = X - 0.5f*limX + 0.5f;
	float fY = Y - 0.5f*limY + 0.5f;
	float fZ = startZ - 0.5f*limZ + 0.5f;

    // solve by small blocks over all angles 
	float Z[ZSIZE];
	for(int i=0; i < ZSIZE; i++)
		Z[i] = 0.0f;


    float fAngle = startAngle + angleOffset + 0.5f;
    float4 fCu, fCv; 
    float fU, fV;
    float fXn, fYn, fZn;  // normalized coordinates 
    float fXs, fYs, fZs;  // shifted coordinates 
    float angle_ratio ; // ratio from angle / iProjAngles 

    for (int angle = startAngle; angle < endAngle; ++angle, fAngle += 1.0f)
    {

        fCu = make_float4(gC_C[8*angle+0], gC_C[8*angle+1], gC_C[8*angle+2], gC_C[8*angle+3]);
        fCv = make_float4(gC_C[8*angle+4], gC_C[8*angle+5], gC_C[8*angle+6], gC_C[8*angle+7]);

        angle_ratio = (float)angle / (float)dims.iProjAngles ; 


        if (use_deform) 
             {

            /*
               // FASTER APPROXIMATION FOR SMALL DEFORMATIONS 
                fXn = X/limX; // normalized coordinates 
                fYn = Y/limY; 
                fZn = startZ/limZ; 
                 

                // load deformed coordinates 
                    fXs =  fX + tex3D(Xdef0_tex,fXn, fYn, fZn); 
                    fYs =  fY + tex3D(Ydef0_tex,fXn, fYn, fZn);
                    fZs =  fZ + tex3D(Zdef0_tex,fXn, fYn, fZn);
                
                // find location on the detector 
                fU = fCu.w + fXs * fCu.x + fYs * fCu.y + fZs * fCu.z;
                fV = fCv.w + fXs * fCv.x + fYs * fCv.y + fZs * fCv.z;

                for (int idx = 0; idx < ZSIZE; ++idx) {
                     // get bilinear interpolation back to non-shifted coordinates 
                    Z[idx] += tex3D(gT_par3DProjTexture, fU, fAngle, fV);

                    // TODO: check if approximation that deformation is constant for Z block is valid  !!
                    fU += fCu.z;
                    fV += fCv.z;
                }
            */

            // ARBITRARY DEFORMATIONS APPROXIMATION

            fXn = X/limX; // normalized coordinates 
            fYn = Y/limY; 
            for (int idx = 0; idx < ZSIZE; ++idx) {
                fZs = fZ + idx; // Z coordinate
                fZn = (startZ+idx)/limZ; //  normalized Z coordinate

                // load deformed coordinates 
                if (!linear_deform_model){ 
                    fXs =  fX + tex3D(Xdef0_tex,fXn, fYn, fZn); 
                    fYs =  fY + tex3D(Ydef0_tex,fXn, fYn, fZn);
                    fZs =  fZs +tex3D(Zdef0_tex,fXn, fYn, fZn);
                } else {
                    // deformated coordinates with linear interpolation 
                    fXs =  fX + (tex3D(Xdef0_tex,fXn, fYn, fZn) * (1-angle_ratio) + angle_ratio*tex3D(Xdef1_tex,fXn, fYn, fZn)); 
                    fYs =  fY + (tex3D(Ydef0_tex,fXn, fYn, fZn) * (1-angle_ratio) + angle_ratio*tex3D(Ydef1_tex,fXn, fYn, fZn));
                    fZs =  fZs +(tex3D(Zdef0_tex,fXn, fYn, fZn) * (1-angle_ratio) + angle_ratio*tex3D(Zdef1_tex,fXn, fYn, fZn));
                }
                
                // find location on the detector 
                fU = fCu.w + fXs * fCu.x + fYs * fCu.y + fZs * fCu.z;
                fV = fCv.w + fXs * fCv.x + fYs * fCv.y + fZs * fCv.z;

                 // get bilinear interpolation back to non-shifted coordinates 
                Z[idx] += tex3D(gT_par3DProjTexture, fU, fAngle, fV);

                // TODO: check if approximation that deformation is constant for Z block is valid  !!
                fU += fCu.z;
                fV += fCv.z;
            }



        } else {


            fU = fCu.w + fX * fCu.x + fY * fCu.y + fZ * fCu.z;
            fV = fCv.w + fX * fCv.x + fY * fCv.y + fZ * fCv.z;

            for (int idx = 0; idx < ZSIZE; ++idx) {

                Z[idx] += tex3D(gT_par3DProjTexture, fU, fAngle, fV);

                fU += fCu.z;
                fV += fCv.z;
            }
        }

    }


	int endZ = ZSIZE;
	if (endZ > dims.iVolZ - startZ)
		endZ = dims.iVolZ - startZ;

	for(int i=0; i < endZ; i++)
		volData[((startZ+i)*dims.iVolY+Y)*volPitch+X] += Z[i] * fOutputScale;
}

// supersampling version
__global__ void dev_par3D_BP_SS(void* D_volData, unsigned int volPitch, int startAngle, int angleOffset, const SDimensions3D dims, float fOutputScale)
{
	float* volData = (float*)D_volData;

	int endAngle = startAngle + g_anglesPerBlock;
	if (endAngle > dims.iProjAngles - angleOffset)
		endAngle = dims.iProjAngles - angleOffset;

	// threadIdx: x = rel x
	//            y = rel y

	// blockIdx:  x = x + y
    //            y = z


	// TO TRY: precompute part of detector intersection formulas in shared mem?
	// TO TRY: inner loop over z, gather ray values in shared mem

	const int X = blockIdx.x % ((dims.iVolX+g_volBlockX-1)/g_volBlockX) * g_volBlockX + threadIdx.x;
	const int Y = blockIdx.x / ((dims.iVolX+g_volBlockX-1)/g_volBlockX) * g_volBlockY + threadIdx.y;

	if (X >= dims.iVolX)
		return;
	if (Y >= dims.iVolY)
		return;

	const int startZ = blockIdx.y * g_volBlockZ;
	int endZ = startZ + g_volBlockZ;
	if (endZ > dims.iVolZ)
		endZ = dims.iVolZ;

	float fX = X - 0.5f*dims.iVolX + 0.5f - 0.5f + 0.5f/dims.iRaysPerVoxelDim;
	float fY = Y - 0.5f*dims.iVolY + 0.5f - 0.5f + 0.5f/dims.iRaysPerVoxelDim;
	float fZ = startZ - 0.5f*dims.iVolZ + 0.5f - 0.5f + 0.5f/dims.iRaysPerVoxelDim;

	const float fSubStep = 1.0f/dims.iRaysPerVoxelDim;

	fOutputScale /= (dims.iRaysPerVoxelDim*dims.iRaysPerVoxelDim*dims.iRaysPerVoxelDim);


	for (int Z = startZ; Z < endZ; ++Z, fZ += 1.0f)
	{

		float fVal = 0.0f;
		float fAngle = startAngle + angleOffset + 0.5f;

		for (int angle = startAngle; angle < endAngle; ++angle, fAngle += 1.0f)
		{
			const float fCux = gC_C[8*angle+0];
			const float fCuy = gC_C[8*angle+1];
			const float fCuz = gC_C[8*angle+2];
			const float fCuc = gC_C[8*angle+3];
			const float fCvx = gC_C[8*angle+4];
			const float fCvy = gC_C[8*angle+5];
			const float fCvz = gC_C[8*angle+6];
			const float fCvc = gC_C[8*angle+7];

			float fXs = fX;
			for (int iSubX = 0; iSubX < dims.iRaysPerVoxelDim; ++iSubX) {
			float fYs = fY;
			for (int iSubY = 0; iSubY < dims.iRaysPerVoxelDim; ++iSubY) {
			float fZs = fZ;
			for (int iSubZ = 0; iSubZ < dims.iRaysPerVoxelDim; ++iSubZ) {

				const float fU = fCuc + fXs * fCux + fYs * fCuy + fZs * fCuz;
				const float fV = fCvc + fXs * fCvx + fYs * fCvy + fZs * fCvz;

				fVal += tex3D(gT_par3DProjTexture, fU, fAngle, fV);
				fZs += fSubStep;
			}
			fYs += fSubStep;
			}
			fXs += fSubStep;
			}

		}

		volData[(Z*dims.iVolY+Y)*volPitch+X] += fVal * fOutputScale;
	}

}

bool Par3DBP_Array(cudaPitchedPtr D_volumeData,
                   const SDimensions3D& dims, const SPar3DProjection* angles,
                   float fOutputScale, bool use_deform, bool linear_deform_model)
{

	for (unsigned int th = 0; th < dims.iProjAngles; th += g_MaxAngles) {
		unsigned int angleCount = g_MaxAngles;
		if (th + angleCount > dims.iProjAngles)
			angleCount = dims.iProjAngles - th;

		// transfer angles to constant memory
		float* tmp = new float[8*dims.iProjAngles];

		// NB: We increment angles at the end of the loop body.


		// TODO: Use functions from dims3d.cu for this:

#define TRANSFER_TO_CONSTANT(expr,name) do { for (unsigned int i = 0; i < angleCount; ++i) tmp[8*i + name] = (expr) ; } while (0)

#define DENOM (angles[i].fRayX*angles[i].fDetUY*angles[i].fDetVZ - angles[i].fRayX*angles[i].fDetUZ*angles[i].fDetVY - angles[i].fRayY*angles[i].fDetUX*angles[i].fDetVZ + angles[i].fRayY*angles[i].fDetUZ*angles[i].fDetVX + angles[i].fRayZ*angles[i].fDetUX*angles[i].fDetVY - angles[i].fRayZ*angles[i].fDetUY*angles[i].fDetVX)

		TRANSFER_TO_CONSTANT( ( - (angles[i].fRayY*angles[i].fDetVZ - angles[i].fRayZ*angles[i].fDetVY)) / DENOM , 0 );
		TRANSFER_TO_CONSTANT( ( (angles[i].fRayX*angles[i].fDetVZ - angles[i].fRayZ*angles[i].fDetVX)) / DENOM , 1 );
		TRANSFER_TO_CONSTANT( (- (angles[i].fRayX*angles[i].fDetVY - angles[i].fRayY*angles[i].fDetVX) ) / DENOM , 2 );
		TRANSFER_TO_CONSTANT( (-(angles[i].fDetSY*angles[i].fDetVZ - angles[i].fDetSZ*angles[i].fDetVY)*angles[i].fRayX + (angles[i].fRayY*angles[i].fDetVZ - angles[i].fRayZ*angles[i].fDetVY)*angles[i].fDetSX - (angles[i].fRayY*angles[i].fDetSZ - angles[i].fRayZ*angles[i].fDetSY)*angles[i].fDetVX) / DENOM , 3 );

		TRANSFER_TO_CONSTANT( ((angles[i].fRayY*angles[i].fDetUZ - angles[i].fRayZ*angles[i].fDetUY) ) / DENOM , 4 );
		TRANSFER_TO_CONSTANT( (- (angles[i].fRayX*angles[i].fDetUZ - angles[i].fRayZ*angles[i].fDetUX) ) / DENOM , 5 );
		TRANSFER_TO_CONSTANT( ((angles[i].fRayX*angles[i].fDetUY - angles[i].fRayY*angles[i].fDetUX) ) / DENOM , 6 );
		TRANSFER_TO_CONSTANT( ((angles[i].fDetSY*angles[i].fDetUZ - angles[i].fDetSZ*angles[i].fDetUY)*angles[i].fRayX - (angles[i].fRayY*angles[i].fDetUZ - angles[i].fRayZ*angles[i].fDetUY)*angles[i].fDetSX + (angles[i].fRayY*angles[i].fDetSZ - angles[i].fRayZ*angles[i].fDetSY)*angles[i].fDetUX ) / DENOM , 7 );

#undef TRANSFER_TO_CONSTANT
#undef DENOM

		cudaMemcpyToSymbol(gC_C, tmp, angleCount*8*sizeof(float), 0, cudaMemcpyHostToDevice); 

		delete[] tmp;

		checkLastError("after cudaMemcpyToSymbol");


		dim3 dimBlock(g_volBlockX, g_volBlockY);

		dim3 dimGrid(((dims.iVolX+g_volBlockX-1)/g_volBlockX)*((dims.iVolY+g_volBlockY-1)/g_volBlockY), (dims.iVolZ+g_volBlockZ-1)/g_volBlockZ);

		// timeval t;
		// tic(t);

		for (unsigned int i = 0; i < angleCount; i += g_anglesPerBlock) {
			// printf("Calling BP: %d, %dx%d, %dx%d to %p\n", i, dimBlock.x, dimBlock.y, dimGrid.x, dimGrid.y, (void*)D_volumeData.ptr); 
			if (dims.iRaysPerVoxelDim == 1)
				dev_par3D_BP<<<dimGrid, dimBlock>>>(D_volumeData.ptr, D_volumeData.pitch/sizeof(float), i, th, dims, fOutputScale, use_deform, linear_deform_model);
			else
				dev_par3D_BP_SS<<<dimGrid, dimBlock>>>(D_volumeData.ptr, D_volumeData.pitch/sizeof(float), i, th, dims, fOutputScale);
		}

		cudaTextForceKernelsCompletion();
		checkLastError("after cudaTextForceKernelsCompletion");

		angles = angles + angleCount;
		// printf("%f\n", toc(t));

	}



	return true;
}

bool Par3DBP(cudaPitchedPtr D_volumeData,
            cudaPitchedPtr D_projData,
            const SDimensions3D& dims, const SPar3DProjection* angles,
            float fOutputScale, DeformField DF)
{
	// transfer projections to array

	checkLastError("before allocateVolumeArray");

	cudaArray* cuArray = allocateProjectionArray(dims);
	checkLastError("after allocateVolumeArray");

	transferProjectionsToArray(D_projData, cuArray, dims);

	checkLastError("after transferProjectionsToArray");

	bindDataTexture(cuArray, gT_par3DProjTexture, cudaAddressModeBorder, false);
	checkLastError("after bindProjDataTexture");

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

	bool ret = Par3DBP_Array(D_volumeData, dims, angles, fOutputScale, DF.use_deform, DF.use_linear_model);

	checkLastError("after Par3DBP_Array");

	cudaUnbindTexture(gT_par3DProjTexture);
	checkLastError("after cudaUnbindTexture");

	cudaFreeArray(cuArray);

	checkLastError("after cudaFreeArray");


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


}

