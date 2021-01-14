/*  
 *
 **-----------------------------------------------------------------------*
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



 *
   Compilation from Matlab:
   mex -R2018a 'CFLAGS="\$CFLAGS -fopenmp"' LDFLAGS="\$LDFLAGS -fopenmp" get_from_3D_projection_mex.cpp
 
  % Usage from Matlab:
  full_array = (randn(1000, 1000, 1, 'single'));
  small_array = (ones(500, 500, 100, 'single'));

  positions = int32([1:100; 1:100])';
  indices = int32([1:100]);  % indices are starting from 1 !! 
  tic; get_from_3D_projection_mex(small_array,full_array,positions,indices); toc

 This code in matlab:
    full_array = randn(100,100,200, 'single'); 
    small_array = zeros(50,50,50, 'single'); 
    positions = ones(200,2, 'int32'); 
    indices = int32(1:50); 
 
    Npix = size(full_array);
    small_array = zeros(dimensions, 'single');
    for jj = 1:length(indices)
        ii = indices(jj)
        for i = 1:2
            % limit to the region inside full_array
            ind_f{i} = max(1,1+positions(ii,i)):min(positions(ii,i)+dimensions(i),Npix(i));
            % adjust size of the small matrix to correspond
            ind_s{i} = ((ind_f{i}(1)-positions(ii,i))):(ind_f{i}(end)-positions(ii,i));
        end
        small_array(ind_s{:},jj) = full_array(ind_f{:},ii) ;
    end
 */

#include "matlab_overload.h"
#include "mex.h"
#include <math.h> 
#include <stdio.h> 
#include <omp.h>
#include <stdint.h>
#include <sys/sysinfo.h>

#define THREADS 12
#define CHUNK 20



template <typename dtype>
void inner_loop(dtype const * array_full, dtype * array_small, const mwSize pos_x0, const mwSize pos_y0, const mwSize pos_zs, const mwSize pos_zf, const mwSize Ns_x, const mwSize Ns_y, const mwSize Nf_x, const mwSize Nf_y)
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

            //skip positions that are out of the matrix
            id_small = row + idc_small;
            id_large = pos_y + idc_large;

            //rewrite original values 
            SetData(array_small[id_small], array_full[id_large]);
         }
    }
}


template <typename dtype>
void get_from_projection(int nlhs, mxArray *plhs[],
                         int nrhs, const mxArray *prhs[])
{
    
  


    dtype const * array_full;  
    dtype       * array_small;

    GetData(prhs[1], array_full);
    GetData(prhs[0], array_small);
     
    // check if values should be added or overwritten 
    bool const add_values = !((nrhs == 5) && ( !mxGetScalar(prhs[4]) ));

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

  if(Npos != Ns_z )
    mexErrMsgIdAndTxt("MexError:tomo","The 3rd dim of 1st input argument must be equal to length of positions.");

  if(Nid != Ns_z) 
    mexErrMsgIdAndTxt("MexError:tomo","The 3rd dim of 1st input argument must be equal to length of indices.");
          
    mwSize id, pos_zs,pos_zf, col, row, pos_y, pos_x, pos_x0, pos_y0, idc_small, idc_large;
    mwSize id_small, id_large;
    bool out_of_range = false;
        
//     #pragma omp parallel for schedule(dynamic) num_threads(THREADS) private(col, row, pos_x, pos_y, pos_x0, pos_y0, id_small, id_large, pos_zs,pos_zf,id, idc_small, idc_large)
    for (id = 0; id < Nid; id++) {
        if (Nf_z == 1)
            pos_zf = 0; 
        else
            pos_zf =  indices[id]-1; // distribute the small_array only to defined sliced in the full_array 
        
        if (pos_zf > Nf_z)
        {
            out_of_range = true; 
            continue;
        }
        pos_zs =  (id < Ns_z ? id : Ns_z);  // min(id, Nf_z) 
        
                
        pos_x0 = positions[id+Nid]; 
        pos_y0 = positions[id]; 

        inner_loop<dtype>(array_full, array_small, pos_x0, pos_y0, pos_zs, pos_zf, Ns_x, Ns_y, Nf_x, Nf_y);

    }
    if (out_of_range)
        mexErrMsgIdAndTxt("MexError:tomo","Indices are out of range for provided inputs");


}


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    #if MX_HAS_INTERLEAVED_COMPLEX == 0
        mexErrMsgIdAndTxt("MexError:tomo","Only Matlab R2018a and newer is supported");
    #endif
   
    /* Check for proper number of arguments. */
    if (nrhs <4 || nrhs > 5)
        mexErrMsgTxt("4-5 input arguments required: add_to_3D_projection_mex(small_array,full_array,positions, indices, add_values)");
    else if (nlhs != 0)
        mexErrMsgTxt("No output argument has to be specified.");
    
    /* Input must be of type double / single / uint32 / uint16. */
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
            case mxDOUBLE_CLASS: get_from_projection<mxComplexDouble>(nlhs, plhs, nrhs, prhs); break; 
            case mxSINGLE_CLASS: get_from_projection<mxComplexSingle>(nlhs, plhs, nrhs, prhs); break; 
            case mxUINT32_CLASS: get_from_projection<mxComplexUint32>(nlhs, plhs, nrhs, prhs); break; 
            case mxUINT16_CLASS: get_from_projection<mxComplexUint16>(nlhs, plhs, nrhs, prhs); break; 
            case mxUINT8_CLASS:  get_from_projection<mxComplexUint8>(nlhs, plhs, nrhs, prhs); break; 
        }
    else
        switch (mxGetClassID(prhs[0]))
        {
            case mxDOUBLE_CLASS: get_from_projection<mxDouble>(nlhs, plhs, nrhs, prhs); break; 
            case mxSINGLE_CLASS: get_from_projection<mxSingle>(nlhs, plhs, nrhs, prhs); break; 
            case mxUINT32_CLASS: get_from_projection<mxUint32>(nlhs, plhs, nrhs, prhs); break; 
            case mxUINT16_CLASS: get_from_projection<mxUint16>(nlhs, plhs, nrhs, prhs); break; 
            case mxUINT8_CLASS:  get_from_projection<mxUint8>(nlhs, plhs, nrhs, prhs); break; 
        }
           
  return;
}
