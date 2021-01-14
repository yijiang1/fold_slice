/*
 *
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


 * Compilation from Matlab:
  mex -R2018a 'CFLAGS="\$CFLAGS -fopenmp"' LDFLAGS="\$LDFLAGS -fopenmp" add_to_3D_projection_mex.cpp
 *
 * Usage from Matlab:
 
  full_array = (zeros(1000, 1000, 1, 'single'));
  small_array = (ones(500, 500, 100, 'single'));

  positions = int32([1:100; 1:100])';
  indices = int32([1:100]);  % indices are starting from 1 !! 
  add_values = true; % (DEFAULT)
  add_to_3D_projection_mex(small_array,full_array,positions, indices,add_values);
 
 * Matlab version: add_to_3D_projection(full_array, small_array, positions)
 *
 *
 *
 * results are directly added to full_array, add_values == false => rewrite original values 
 *
 * This code in matlab:
 *
  N_f = size(full_array);
  N_s = size(small_array);
  for ii = 1:N_f(3)
  for i = 1:2
  ind_f{i} = max(1, 1+positions(ii,i)):min(N_f(i),positions(ii,i)+N_s(i));
  ind_s{i} = ((ind_f{i}(1)-positions(ii,i))):(ind_f{i}(end)-positions(ii,i));
  end
  full_array(ind_f{:},ii) = full_array(ind_f{:},ii) + small_array(ind_s{:},ii);
  end
 *
 */

#include "matlab_overload.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <omp.h>
#include <stdint.h>


#define THREADS 16
#define CHUNK 10

template <typename dtype, bool add_atomic, bool add_values>
void inner_loop(dtype * array_full, dtype const * array_small, const mwSize pos_x0, const mwSize pos_y0, const mwSize pos_zs, const mwSize pos_zf, const mwSize Ns_x, const mwSize Ns_y, const mwSize Nf_x, const mwSize Nf_y)
{
    
    mwSize id, col, row, pos_y, pos_x, id_small, id_large, idc_small, idc_large;
    
    #pragma omp parallel for schedule(static) num_threads(THREADS) private(col, row, pos_x, pos_y, id_small, id_large, idc_small, idc_large)
    for (col = (pos_x0>=0 ? 0 : -pos_x0) ; col < Ns_x; col++) {
        pos_x = col + pos_x0;
        idc_small = col*Ns_y + Ns_y*Ns_x*pos_zs; 
        idc_large = pos_x*Nf_y + Nf_y*Nf_x*pos_zf; 

        if (pos_x >= Nf_x )
            continue; 
        for (row = (pos_y0 >= 0 ? 0 : -pos_y0) ; row < Ns_y; row++) {
            pos_y = row + pos_y0;
            if (pos_y >= Nf_y )
                continue; 

            // skip positions that are out of the matrix
            id_small = row + idc_small;
            id_large = pos_y + idc_large;

            if (add_atomic && add_values)
                //Add values to the already provided ones 
                 AddData_atomic(array_full[id_large], array_small[id_small]);
            else if (add_values) 
                // rewrite original values 
                 AddData(array_full[id_large], array_small[id_small]);
            else
                // rewrite original values 
                 SetData(array_full[id_large], array_small[id_small]);
        }
    }
}


template <typename dtype>
void add_to_projection(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{


    dtype       * array_full;  
    dtype const * array_small;

    GetData(prhs[1], array_full);
    GetData(prhs[0], array_small);
     
    // check if values should be added or overwritten 
    bool const add_values = nrhs < 5 || mxGetScalar(prhs[4]);  // if true, x += y, if false x = y; 
    bool const add_atomic = nrhs < 6 || mxGetScalar(prhs[5]);  // if true, correclty deal with overlap between the positions, but it is slow 

    mxInt32 const *indices =   mxGetInt32s(prhs[3]);
    mxInt32 const *positions = mxGetInt32s(prhs[2]);
    
    /* Get dimension of probe and object / small + large array  */
    mwSize const *fdims = mxGetDimensions(prhs[1]);
    mwSize const Nf_y = fdims[0];
    mwSize const Nf_x = fdims[1];
    mwSize const Nf_z = (mxGetNumberOfDimensions(prhs[1]) == 3 ? fdims[2] : 1); 
    
    mwSize const *sdims = mxGetDimensions(prhs[0]);
    mwSize const Ns_y = sdims[0];
    mwSize const Ns_x = sdims[1];
    mwSize const Ns_z = (mxGetNumberOfDimensions(prhs[0]) == 3 ? sdims[2] : 1); 
     
    mwSize const Nid = mxGetNumberOfElements(prhs[3]);
    mwSize const Npos = mxGetM(prhs[2]);
 
  if(Npos != Ns_z && Ns_z > 1 )
    mexErrMsgIdAndTxt("MexError:tomo","The 3rd dim of 1st input argument must be equal to length of positions.");

  if(Nid != Ns_z && Ns_z > 1) 
    mexErrMsgIdAndTxt("MexError:tomo","The 3rd dim of 1st input argument must be equal to length of indices.");
          
    mwSize id, pos_zs,pos_zf, col, row, pos_y, pos_x, pos_x0, pos_y0 ; 
    mwSize id_small, id_large, idc_small, idc_large;
    bool out_of_range = false; 
    
    for (id = 0; id < Nid; id++) {
        if (Nf_z == 1)
            pos_zf = 0; 
        else
            pos_zf =  indices[id]-1; // distribute the small_array only to defined sliced in the full_array 
                
        if (pos_zf >= Nf_z)
        {
            out_of_range = true; 
            continue;
        }
        pos_zs = (id < Ns_z ? id : Ns_z-1);  // min(id, Nf_z) 
        
                
        pos_x0 = positions[id+Nid]; 
        pos_y0 = positions[id]; 
       
        if (add_values && add_atomic)
            inner_loop<dtype,true,true>(array_full, array_small, pos_x0, pos_y0, pos_zs, pos_zf, Ns_x, Ns_y, Nf_x, Nf_y);
        else if (add_values)
            inner_loop<dtype,false,true>(array_full, array_small, pos_x0, pos_y0, pos_zs, pos_zf, Ns_x, Ns_y, Nf_x, Nf_y);
        else
            inner_loop<dtype,false,false>(array_full, array_small, pos_x0, pos_y0, pos_zs, pos_zf, Ns_x, Ns_y, Nf_x, Nf_y);

    }
    if (out_of_range)
        mexErrMsgIdAndTxt("MexError:tomo","Indices are out of range for provided inputs");

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    #if MX_HAS_INTERLEAVED_COMPLEX == 0
        mexErrMsgIdAndTxt("MexError:tomo","Only Matlab R2018a and newer is supported");
    #endif
   
    /* Check for proper number of arguments. */
    if (nrhs <4 || nrhs > 6)
        mexErrMsgTxt("4-6 input arguments required: add_to_3D_projection_mex(small_array,full_array,positions, indices, add_values=true, add_atomic=true)");
    else if (nlhs != 0)
        mexErrMsgTxt("No output argument has to be specified.");
    
    /* Input must be of type single / uint32 / uint16. */
    if ( !(mxIsDouble(prhs[0]) || mxIsSingle(prhs[0]) || mxIsUint32(prhs[0]) || mxIsUint16(prhs[0]) || mxIsLogical(prhs[0]) || mxIsUint8(prhs[0]) ) ) {
        mexErrMsgIdAndTxt("MexError:tomo","Class of input 1 is not double/single/uint8/uint16/uint32");
    }
    if (  (mxGetClassID(prhs[0]) != mxGetClassID (prhs[1])) ) {
        mexErrMsgIdAndTxt("MexError:tomo","Inputs arrays are not the same type");
    }
    /* Input must be of type int32. */
    for (int i=2; i<4; i++) {
        if (mxIsInt32(prhs[i]) != 1) {
            printf("Input %d is not integer\n",i+1);
            mexErrMsgIdAndTxt("MexError:tomo","Inputs must be of correct type.");
        }
    }
    if ((nrhs == 5) && (mxIsLogical(prhs[4]) != 1)) {
        printf("Input 5 is not logical\n");
        mexErrMsgIdAndTxt("MexError:tomo","Inputs must be of correct type.");
    }
    
    
    if(mxIsComplex(prhs[0]) !=  mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MexError:tomo","Complexity of the inputs has to be the same");
    }
    
    
    
  if ((mxGetNumberOfDimensions(prhs[0]) > 3) || (mxGetNumberOfDimensions(prhs[0]) < 2) ||
      (mxGetNumberOfDimensions(prhs[1]) > 3) || (mxGetNumberOfDimensions(prhs[1]) < 2) ||
      (mxGetNumberOfDimensions(prhs[2]) != 2) ||
      (mxGetNumberOfDimensions(prhs[3]) != 2))
      mexErrMsgIdAndTxt("MexError:tomo","Wrong number of dimensions in inputs");
    
    
  if(mxGetN(prhs[2]) != 2 )
    mexErrMsgIdAndTxt("MexError:tomo","Positions are expected as Nx2 matrix");
        
   if (mxIsComplex(prhs[0]))
        switch (mxGetClassID(prhs[0]))
        {
            case mxDOUBLE_CLASS: add_to_projection<mxComplexDouble>(nlhs, plhs, nrhs, prhs); break; 
            case mxSINGLE_CLASS: add_to_projection<mxComplexSingle>(nlhs, plhs, nrhs, prhs); break; 
            case mxUINT32_CLASS: add_to_projection<mxComplexUint32>(nlhs, plhs, nrhs, prhs); break; 
            case mxUINT16_CLASS: add_to_projection<mxComplexUint16>(nlhs, plhs, nrhs, prhs); break; 
            case mxUINT8_CLASS:  add_to_projection<mxComplexUint8>(nlhs, plhs, nrhs, prhs); break; 
        }
    else
        switch (mxGetClassID(prhs[0]))
        {
            case mxDOUBLE_CLASS: add_to_projection<mxDouble>(nlhs, plhs, nrhs, prhs); break; 
            case mxSINGLE_CLASS: add_to_projection<mxSingle>(nlhs, plhs, nrhs, prhs); break; 
            case mxUINT32_CLASS: add_to_projection<mxUint32>(nlhs, plhs, nrhs, prhs); break; 
            case mxUINT16_CLASS: add_to_projection<mxUint16>(nlhs, plhs, nrhs, prhs); break; 
            case mxUINT8_CLASS:  add_to_projection<mxUint8>(nlhs, plhs, nrhs, prhs); break; 
        }
    
}

