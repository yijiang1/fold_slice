/*  
 *
 *Academic License Agreement

Source Code

Introduction 
•	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
  will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the cSAXS 
  ptychography MATLAB package computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").

Terms and Conditions of the LICENSE
1.	LICENSOR grants to LICENSEE a royalty-free, non-exclusive license to use the PROGRAM for academic, non-commercial purposes, upon the terms and conditions 
      hereinafter set out and until termination of this license as set forth below.
2.	LICENSEE acknowledges that the PROGRAM is a research tool still in the development stage. The PROGRAM is provided without any related services, improvements 
      or warranties from LICENSOR and that the LICENSE is entered into in order to enable others to utilize the PROGRAM in their academic activities. It is the 
      LICENSEE’s responsibility to ensure its proper use and the correctness of the results.”
3.	THE PROGRAM IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR 
      A PARTICULAR PURPOSE AND NONINFRINGEMENT OF ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. IN NO EVENT SHALL THE LICENSOR, THE AUTHORS OR THE COPYRIGHT 
      HOLDERS BE LIABLE FOR ANY CLAIM, DIRECT, INDIRECT OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION WITH THE PROGRAM OR THE USE 
      OF THE PROGRAM OR OTHER DEALINGS IN THE PROGRAM.
4.	LICENSEE agrees that it will use the PROGRAM and any modifications, improvements, or derivatives of PROGRAM that LICENSEE may create (collectively, 
      "IMPROVEMENTS") solely for academic, non-commercial purposes and that any copy of PROGRAM or derivatives thereof shall be distributed only under the same 
      license as PROGRAM. The terms "academic, non-commercial", as used in this Agreement, mean academic or other scholarly research which (a) is not undertaken for 
      profit, or (b) is not intended to produce works, services, or data for commercial use, or (c) is neither conducted, nor funded, by a person or an entity engaged 
      in the commercial use, application or exploitation of works similar to the PROGRAM.
5.	LICENSEE agrees that it shall make the following acknowledgement in any publication resulting from the use of the PROGRAM or any translation of the code into 
      another computing language:
      "Data processing was carried out using the cSAXS ptychography MATLAB package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
      Scherrer Institut, Switzerland."

Additionally, any publication using the package, or any translation of the code into another computing language should cite for difference map:
P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379–382 (2008). 
  (doi: 10.1126/science.1158573),
for mixed coherent modes:
P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
for LSQ-ML method 
M. Odstrcil, A. Menzel, M.G. Sicairos,  Iterative least-squares solver for generalized maximum-likelihood ptychography, Optics Express, 2018
for OPRP method 
 M. Odstrcil, P. Baksh, S. A. Boden, R. Card, J. E. Chad, J. G. Frey, W. S. Brocklesby,  "Ptychographic coherent diffractive imaging with orthogonal probe relaxation." Optics express 24.8 (2016): 8360-8369
and/or for multislice:
E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
6.	Except for the above-mentioned acknowledgment, LICENSEE shall not use the PROGRAM title or the names or logos of LICENSOR, nor any adaptation thereof, nor the 
      names of any of its employees or laboratories, in any advertising, promotional or sales material without prior written consent obtained from LICENSOR in each case.
7.	Ownership of all rights, including copyright in the PROGRAM and in any material associated therewith, shall at all times remain with LICENSOR, and LICENSEE 
      agrees to preserve same. LICENSEE agrees not to use any portion of the PROGRAM or of any IMPROVEMENTS in any machine-readable form outside the PROGRAM, nor to 
      make any copies except for its internal use, without prior written consent of LICENSOR. LICENSEE agrees to place the following copyright notice on any such copies: 
      © All rights reserved. PAUL SCHERRER INSTITUT, Switzerland, Laboratory for Macromolecules and Bioimaging, 2017. 
8.	The LICENSE shall not be construed to confer any rights upon LICENSEE by implication or otherwise except as specifically set forth herein.
9.	DISCLAIMER: LICENSEE shall be aware that Phase Focus Limited of Sheffield, UK has an international portfolio of patents and pending applications which relate 
      to ptychography and that the PROGRAM may be capable of being used in circumstances which may fall within the claims of one or more of the Phase Focus patents, 
      in particular of patent with international application number PCT/GB2005/001464. The LICENSOR explicitly declares not to indemnify the users of the software 
      in case Phase Focus or any other third party will open a legal action against the LICENSEE due to the use of the program.
10.	This Agreement shall be governed by the material laws of Switzerland and any dispute arising out of this Agreement or use of the PROGRAM shall be brought before 
      the courts of Zürich, Switzerland. 


 *
 *
   Compilation from Matlab:
   maybe a tiny bit faster code is generated by 
   mex -O COPTIMFLAGS='-O2' LDOPTIMFLAGS='-O2' set_views_cpu_mex.c
 
   Usage from Matlab:
   set_views_cpu_mex(probe,object,positions, Npos);
 
 This code in matlab:
 asize = size(probe);
 for i=1:Npos
    Indy = positions(i,1) + (1:asize(1));
    Indx = positions(i,2) + (1:asize(2));
    ob(Indy,Indx) = ob(Indy,Indx) + probe;
 end
 */

#include "mex.h"
#include <math.h> 
#include <stdio.h> 
#include <omp.h>


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  int i;
    
    /* Check for proper number of arguments. */
  if (nrhs != 4) 
    mexErrMsgTxt("Four input arguments required: set_views_cpu_mex(probe,object,positions,Npos)");
  else if (nlhs != 0) 
    mexErrMsgTxt("No output argument has to be specified.");
  
    /* Input must be of type single. */
  for (i=0; i < 2; i++) {
      if (mxIsSingle(prhs[i]) != 1){
          printf("Input %d is not single\n",i+1);
           mexErrMsgIdAndTxt("MexError:ptycho","Inputs must be of correct type.");
      }
  }
    /* Input must be of type int32. */
  for (i=2; i<nrhs; i++){
    if (mxIsInt32(prhs[i]) != 1){
          printf("Input %d is not integer\n",i+1);
          mexErrMsgIdAndTxt("MexError:ptycho","Inputs must be of correct type.");
    }
  }
 
    /* It cannot be one-dimensional */
  if(mxGetNumberOfDimensions(prhs[0]) < 2) {
    printf("The 1st input argument must have at least two dimensions.");
    mexErrMsgIdAndTxt("MexError:ptycho","wrong number of dimensions");
  }
    /* It cannot be more than 3-dimensional */
   if(mxGetNumberOfDimensions(prhs[0]) > 3) {
    printf("The 1st input argument must have at most three dimensions.");
    mexErrMsgIdAndTxt("MexError:ptycho","wrong number of dimensions");
  }
//   /* Check that arrays are complex */
//   if(mxIsComplex(prhs[0]) != 1) {
//       printf("object input argument must be complex-valued.");
//       mexErrMsgIdAndTxt("MexError:ptycho","Expected complex arrays");
//   }
//   if(mxIsComplex(prhs[1]) != 1) {
//       printf("probe input argument must be complex-valued.");
//       mexErrMsgIdAndTxt("MexError:ptycho","Expected complex arrays");
//   } 
  
  float *object_r, *object_i, *probe_r,*probe_i;
  const int *positions, *ind_ok;    
  bool cprobe, cobject;
          
  cobject =  mxIsComplex(prhs[0]);
  cprobe  =  mxIsComplex(prhs[1]);
  if( cobject != cobject) 
  {
          printf("probe/object input argument must be complex/real-valued.");
          mexErrMsgIdAndTxt("MexError:ptycho","Expected both complex / real arrays");
  }
 
    ind_ok = (int*)mxGetData(prhs[3]);
    positions = (int*)mxGetData(prhs[2]);
    object_r = (float*)mxGetData(prhs[0]);
    probe_r = (float*)mxGetData(prhs[1]);
    if(cprobe)
    {
        /* get pointers to input data */
        object_i = (float*)mxGetImagData(prhs[0]);
        probe_i = (float*)mxGetImagData(prhs[1]);
    }

  
  /* Get dimension of probe and object */
  const mwSize Ndims = mxGetNumberOfDimensions(prhs[1]);
  const mwSize * dims = mxGetDimensions(prhs[1]);
  const mwSize No_y = mxGetM(prhs[0]);
  const mwSize No_x = mxGetN(prhs[0]);
  const mwSize Np_y = dims[0];
  const mwSize Np_x = dims[1];
  const mwSize Npos = mxGetM(prhs[2]);
  
  if((mxGetM(prhs[2]) != dims[2]) && (Ndims == 3)) {
      printf("wrong size of update / positions  %i", Ndims);
      mexErrMsgIdAndTxt("MexError:ptycho","wrong size of update / positions");
  } 
  
  mwSize id_small, id_large, pos, col, row, o, p; 
  bool flat_probe = Ndims == 2;


  #pragma omp parallel for private(p,pos, col, row, id_small, id_large)
  for (p=0;p<Npos;p++){
      pos = ind_ok[p]-1;
      for (col=0;col<Np_x;col++) {
          for (row=0;row<Np_y;row++) {
              if(flat_probe)
                  id_small = row + col*Np_y;
              else
                  id_small = row + col*Np_y + Np_y*Np_x*pos;

              id_large = row + positions[pos] + (col+positions[pos+Npos])*No_y;
              
              #pragma omp atomic
              object_r[id_large] +=  probe_r[id_small];
              if(cprobe)
                #pragma omp atomic
                object_i[id_large] +=  probe_i[id_small];
          }
      }
  }

  
  return;
}
