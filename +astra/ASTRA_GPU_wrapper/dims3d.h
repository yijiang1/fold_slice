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

#ifndef _CUDA_CONE_DIMS_H
#define _CUDA_CONE_DIMS_H

#include "astra/GeometryUtil3D.h"
#include "mex.h"
#include "gpu/mxGPUArray.h"


namespace astraCUDA3d {

using astra::SConeProjection;
using astra::SPar3DProjection;

struct SDimensions3D {
	unsigned int iVolX;
	unsigned int iVolY;
	unsigned int iVolZ;
	unsigned int iProjAngles;
	unsigned int iProjU; // number of detectors in the U direction
	unsigned int iProjV; // number of detectors in the V direction
	unsigned int iRaysPerDetDim;
	unsigned int iRaysPerVoxelDim;
};

struct DeformField {
	const mxGPUArray * X0;
	const mxGPUArray * Y0;
	const mxGPUArray * Z0;
    const mxGPUArray * X1;
	const mxGPUArray * Y1;
	const mxGPUArray * Z1;
    bool  use_deform;
    bool  use_linear_model;
    
};

}

#endif

