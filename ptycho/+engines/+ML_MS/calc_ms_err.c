/* Compile from matlab with:
mex -largeArrayDims 'CFLAGS="\$CFLAGS -std=c99 -fopenmp"' LDFLAGS="\$LDFLAGS -fopenmp" -lmkl_core -lmkl_intel_ilp64 -lmkl_sequential -lmkl_avx2 calc_ms_err.c

Academic License Agreement

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
 for maximum likelihood:
 P. Thibault and M. Guizar-Sicairos, Maximum-likelihood refinement for coherent diffractive imaging, New J. Phys. 14, 063004 (2012). 
   (doi: 10.1088/1367-2630/14/6/063004),
 for mixed coherent modes:
 P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
 and/or for multislice:
 E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
   (doi: 10.1364/OE.24.029089).
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
*/

#include "mex.h"
#include <mkl.h>
#include <omp.h>

void mexFunction(int n_out, mxArray *p_out[], int n_in, const mxArray *p_in[])
{
    // read input values
    MKL_Complex16 *xopt = (MKL_Complex16 *) mxGetData(p_in[0]);
    
    size_t num_objs = (size_t) mxGetScalar(p_in[1]); // can be either 1 or 2
    size_t num_slices = (size_t) mxGetScalar(p_in[2]);
    
    size_t obj_rows = (size_t) mxGetScalar(p_in[3]);
    size_t obj_cols = (size_t) mxGetScalar(p_in[4]);
    
    size_t scan_num = mxGetM(p_in[5]);
    size_t *scan_pos_row = (size_t *) mxGetData(p_in[5]);
    size_t *scan_pos_col = (size_t *) mxGetData(p_in[6]);
    
    size_t scan_rows = (size_t) mxGetScalar(p_in[7]);
    size_t scan_cols = (size_t) mxGetScalar(p_in[8]);
    
    MKL_Complex16 *prop = (MKL_Complex16 *) mxGetData(p_in[9]);
    
    double fnorm_inv = mxGetScalar(p_in[10]);
    double dz = mxGetScalar(p_in[11]);
    
    char *fmask = (char *) mxGetData(p_in[12]);
    double *fmag = mxGetPr(p_in[13]);
    
    // Decompose xopt into object slices and a probe
    size_t i, j;
    MKL_Complex16 *obj[num_objs][num_slices];
    for (i = 0; i < num_objs; i++) {
        for (j = 0; j < num_slices; j++) {
            obj[i][j] = xopt + (i*num_slices+j)*obj_rows*obj_cols;
        }
    }
    
    MKL_Complex16 *probe = xopt+(num_objs*num_slices)*obj_rows*obj_cols;
    
    // arrange FFT descriptor
    DFTI_DESCRIPTOR_HANDLE desc_handle;
    MKL_LONG l[2];
    l[0] = scan_rows; l[1] = scan_cols;
    DftiCreateDescriptor(&desc_handle, DFTI_DOUBLE, DFTI_COMPLEX, 2, l);
    DftiSetValue(desc_handle, DFTI_COMPLEX_STORAGE, DFTI_COMPLEX_COMPLEX);
    DftiSetValue(desc_handle, DFTI_BACKWARD_SCALE, 1./(scan_rows*scan_cols));
    DftiCommitDescriptor(desc_handle);
    
    // calculate error
    double err = 0.0;
    size_t ind_obj, ind_scan;
    #pragma omp parallel private(j, ind_obj, ind_scan) reduction(+:err)
    {
        MKL_Complex16 *obj_s0, *obj_s1;
        MKL_Complex16 *temp = (MKL_Complex16 *) mkl_malloc(scan_rows*scan_cols*sizeof(MKL_Complex16), 64);
        double *buff = (double *) mkl_malloc(scan_rows*scan_cols*sizeof(double), 64);
        
        // loop over beam positions
        #pragma omp for
        for (i = 0; i < scan_num; i++) {
            ind_obj = scan_pos_col[i]*obj_rows + scan_pos_row[i];
            ind_scan = i*scan_rows*scan_cols;
            
            if (num_objs == 1) {
                obj_s0 = obj[0][0];
                obj_s1 = obj[0][1];
                
            } else {
                obj_s0 = (i < scan_num/2) ? obj[0][0] : obj[1][0];
                obj_s1 = (i < scan_num/2) ? obj[0][1] : obj[1][1];
            }
            
            for (j = 0; j < scan_cols; j++) {
                vzMul(scan_rows, &obj_s0[ind_obj+j*obj_rows], &probe[j*scan_rows], &temp[j*scan_rows]);
            }
            
            if (dz != 0) {
                DftiComputeForward(desc_handle, temp);
                vzMul(scan_rows*scan_cols, temp, prop, temp);
                DftiComputeBackward(desc_handle, temp);
            }
            
            for (j = 0; j < scan_cols; j++) {
                vzMul(scan_rows, &temp[j*scan_rows], &obj_s1[ind_obj+j*obj_rows], &temp[j*scan_rows]);
            }
            
            DftiComputeForward(desc_handle, temp);
            
            vzAbs(scan_rows*scan_cols, temp, buff);
            
            for (j = 0; j < scan_rows*scan_cols; j++) {
                buff[j] = buff[j]*fnorm_inv-fmag[ind_scan+j];
                err += fmask[ind_scan+j]*buff[j]*buff[j];
            }
        }
        
        mkl_free(temp);
        mkl_free(buff);
    }
    
    DftiFreeDescriptor(&desc_handle);
    mkl_free_buffers();
    
    // create output stack of slices
    p_out[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
    double *out = mxGetPr(p_out[0]);
    *out = err;
}

