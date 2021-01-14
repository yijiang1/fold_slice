// RADIAL_INTEG_MEX perform radial integration for a 2D frame
// 
// ** ind_r_max         int32
// ** no_of_segments    int32
// ** norm_sum          double
// ** indices           cell
// ** frame_data        2D or 3D array, double
// 
// return:
// ++ frame_I           2D or 3D array, double
// ++ frame_std         2D or 3D array, double
// 
// 
// Example:
// [frame_I(:,:,ind_frame), frame_std(:,:,ind_frame)] = radial_integ_mex(int32(ind_r_max),int32(no_of_segments), integ_masks.norm_sum, integ_masks.indices, frame_data);
// 
// compile with:
// mex 'CFLAGS="\$CFLAGS -O3 -std=c++17 -fopenmp"' LDFLAGS="\$LDFLAGS -fopenmp" radial_integ_mex.cpp
// 
// MATLAB code:
//    for (ind_r = 1:ind_r_max)
//         for (ind_seg = 1:no_of_segments)
//             if (integ_masks.norm_sum(ind_r,ind_seg) > 0)
//                 frame_I_one_frame(ind_r,ind_seg) = ...
//                     mean(frame_data(integ_masks.indices{ind_r,ind_seg}));
//                 frame_std_one_frame(ind_r,ind_seg) = ...
//                     std(frame_data(integ_masks.indices{ind_r,ind_seg}));
//             else
//                 % mark unknown intensities
//                 frame_I_one_frame(ind_r,ind_seg) = -1;
//                 frame_std_one_frame(ind_r,ind_seg) = -1;
//             end
//         end
//     end

// *-----------------------------------------------------------------------*
// |                                                                       |
// |  Except where otherwise noted, this work is licensed under a          |
// |  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
// |  International (CC BY-NC-SA 4.0) license.                             |
// |                                                                       |
// |  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
// |                                                                       |
// |       Author: CXS group, PSI                                          |
// *-----------------------------------------------------------------------*
// You may use this code with the following provisions:
// 
// If the code is fully or partially redistributed, or rewritten in another
//   computing language this notice should be included in the redistribution.
// 
// If this code, or subfunctions or parts of it, is used for research in a 
//   publication or if it is fully or partially rewritten for another 
//   computing language the authors and institution should be acknowledged 
//   in written form in the publication: “Data processing was carried out 
//   using the “cSAXS matlab package” developed by the CXS group,
//   Paul Scherrer Institut, Switzerland.” 
//   Variations on the latter text can be incorporated upon discussion with 
//   the CXS group if needed to more specifically reflect the use of the package 
//   for the published work.
// 
// A publication that focuses on describing features, or parameters, that
//    are already existing in the code should be first discussed with the
//    authors.
//   
// This code and subroutines are part of a continuous development, they 
//    are provided “as they are” without guarantees or liability on part
//    of PSI or the authors. It is the user responsibility to ensure its 
//    proper use and the correctness of the results.

#include "mex.h"
#include "matrix.h"
#include <iostream>
#include <math.h>
#include <omp.h>

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double *norm_sum;
    double *frame_data;
    uint ind_r_max, no_of_segments;     
    const mwSize *pDims;
    int nDimNum;
    int maxSlice;
    
    /* check for proper number of arguments */
    if(nrhs!=5) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Five inputs required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Two output containers are required.");
    }
    /* make sure the first two input arguments are of type int */
    if( !mxIsClass(prhs[0], "int32")) {
        mexErrMsgIdAndTxt("cxsSoftware:radialIntegMex:notInteger","r_max must be of type integer.");
    }
    if( !mxIsClass(prhs[1], "int32")) {
        mexErrMsgIdAndTxt("cxsSoftware:radialIntegMex:notInteger","no_of_segments must be of type integer.");
    }
    if( !mxIsClass(prhs[2], "double")) {
        mexErrMsgIdAndTxt("cxsSoftware:radialIntegMex:notInteger","norm_sum must be of type double.");
    }
    if( !mxIsCell(prhs[3])) {
        mexErrMsgIdAndTxt("cxsSoftware:radialIntegMex:notInteger","indices must be of type cell.");
    }
    if( !mxIsClass(prhs[4], "double")) {
        mexErrMsgIdAndTxt("cxsSoftware:radialIntegMex:notInteger","frame_data must be of type double.");
    }
    ind_r_max = mxGetScalar(prhs[0]);
    no_of_segments = mxGetScalar(prhs[1]);
    norm_sum = mxGetPr(prhs[2]);
    
    frame_data = mxGetPr(prhs[4]);
    
    nDimNum = mxGetNumberOfDimensions(prhs[4]);
    pDims = mxGetDimensions(prhs[4]);

    if (nDimNum==2){
        plhs[0] = mxCreateNumericMatrix((mwSize)ind_r_max, (mwSize)no_of_segments, mxDOUBLE_CLASS, mxREAL);
        plhs[1] = mxCreateNumericMatrix((mwSize)ind_r_max, (mwSize)no_of_segments, mxDOUBLE_CLASS, mxREAL);
        maxSlice = 1;
    } else if (nDimNum==3) {
        maxSlice = pDims[2];
        mwSize dims[3] = {(mwSize)ind_r_max,(mwSize)no_of_segments,(mwSize)maxSlice};
        plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
        plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    } else {
        mexErrMsgIdAndTxt("cxsSoftware:radialIntegMex:dimsError","frame_data must be 2D or 3D.");
    }
    
    double* outputMatrixMean = (double *)mxGetData(plhs[0]);
    double* outputMatrixStdDev = (double *)mxGetData(plhs[1]);
    
    #pragma omp parallel for collapse(2)
    for (uint slID=0; slID < maxSlice; slID++){
        for (uint ind_r=0; ind_r<ind_r_max; ind_r++){
            for (uint ind_seg=0; ind_seg<no_of_segments; ind_seg++){
                double tmpMean = 0;
                double tmpStdDev = 0;
                
                if (norm_sum[ind_r+ind_seg*ind_r_max] > 0){
                    mxArray *subarray = mxGetCell(prhs[3], ind_r+ind_seg*ind_r_max);
                    double *indicesSub = mxGetPr(subarray);
                    uint dim = mxGetNumberOfElements(subarray);
                    
                    
                    for (uint ii=0; ii<dim; ii++){
                        tmpMean += frame_data[(int)indicesSub[ii]-1 + slID*pDims[0]*pDims[1]];                        
                    }
                    if (tmpMean>0){
                        tmpMean /= dim;
                    }
                    for (uint ii=0; ii<dim; ii++){
                        tmpStdDev += std::pow(std::abs(frame_data[(int)indicesSub[ii]-1 + slID*pDims[0]*pDims[1]]-tmpMean),2);
                        
                    }
                    
                    if (tmpMean>0 && dim>1){
                        tmpStdDev = sqrt(tmpStdDev/(dim-1));
                    }
                    
                    
                } else {
                    tmpMean = -1;
                    tmpStdDev = -1;
                }
                
                outputMatrixMean[ind_r+ind_seg*ind_r_max+slID*(ind_r_max)*(no_of_segments)] = tmpMean;
                outputMatrixStdDev[ind_r+ind_seg*ind_r_max+slID*(ind_r_max)*(no_of_segments)] = tmpStdDev;
            }
        }

    }
        
        


    
}