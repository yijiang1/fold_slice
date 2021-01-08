/******************************************************************************
 *
 * MATLAB (R) is a trademark of The Mathworks (R) Corporation
 *
 * Function:    halfprecision
 * Filename:    halfprecision.c
 * Programmer:  James Tursa
 * Version:     1.0
 * Date:        March 3, 2009
 * Copyright:   (c) 2009 by James Tursa, All Rights Reserved
 *
 *  This code uses the BSD License:
 *
 *  Redistribution and use in source and binary forms, with or without 
 *  modification, are permitted provided that the following conditions are 
 *  met:
 *
 *     * Redistributions of source code must retain the above copyright 
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright 
 *       notice, this list of conditions and the following disclaimer in 
 *       the documentation and/or other materials provided with the distribution
 *      
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * halfprecision converts the input argument to/from a half precision floating
 * point bit pattern corresponding to IEEE 754r. The bit pattern is stored in a
 * uint16 class variable. Please note that halfprecision is *not* a class. That
 * is, you cannot do any arithmetic with the half precision bit patterns.
 * halfprecision is simply a function that converts the IEEE 754r half precision
 * bit pattern to/from other numeric MATLAB variables. You can, however, take
 * the half precision bit patterns, convert them to single or double, do the
 * operation, and then convert the result back manually.
 *
 * 1 bit sign bit
 * 5 bits exponent, biased by 15
 * 10 bits mantissa, hidden leading bit, normalized to 1.0
 *
 * Special floating point bit patterns recognized and supported:
 *
 * All exponent bits zero:
 * - If all mantissa bits are zero, then number is zero (possibly signed)
 * - Otherwise, number is a denormalized bit pattern
 *
 * All exponent bits set to 1:
 * - If all mantissa bits are zero, then number is +Infinity or -Infinity
 * - Otherwise, number is NaN (Not a Number)
 *
 * Building:
 *
 *  halfprecision requires that a mex routine be built (one time only). This
 *  process is typically self-building the first time you call the function
 *  as long as you have the files halfprecision.m and halfprecision.c in the
 *  same directory somewhere on the MATLAB path. If you need to manually build
 *  the mex function, here are the commands:
 *
 * >> mex -setup
 *   (then follow instructions to select a C / C++ compiler of your choice)
 * >> mex halfprecision.c
 *
 * If you have an older version of MATLAB, you may need to use this command:
 *
 * >> mex -DDEFINEMWSIZE halfprecision.c
 *
 * Syntax
 *
 * B = halfprecision(A)
 * C = halfprecision(B,S)
 *     halfprecision(B,'disp')
 *
 * Description
 *
 * A = a MATLAB numeric array, char array, or logical array.
 *
 * B = the variable A converted into half precision floating point bit pattern.
 *     The bit pattern will be returned as a uint16 class variable. The values
 *     displayed are simply the bit pattern interpreted as if it were an unsigned
 *     16-bit integer. To see the halfprecision values, use the 'disp' option, which
 *     simply converts the bit patterns into a single class and then displays them.
 *
 * C = the half precision floating point bit pattern in B converted into class S.
 *     B must be a uint16 or int16 class variable.
 *
 * S = char string naming the desired class (e.g., 'single', 'int32', etc.)
 *     If S = 'disp', then the floating point bit values are simply displayed.
 *
 * Examples
 * 
 * >> a = [-inf -1e30 -1.2 NaN 1.2 1e30 inf]
 * a =
 * 1.0e+030 *
 *     -Inf   -1.0000   -0.0000       NaN    0.0000    1.0000       Inf
 *
 * >> b = halfprecision(a)
 * b =
 * 64512  64512  48333  65024  15565  31744  31744
 *
 * >> halfprecision(b,'disp')
 *     -Inf      -Inf   -1.2002       NaN    1.2002       Inf       Inf
 *
 * >> halfprecision(b,'double')
 * ans =
 *     -Inf      -Inf   -1.2002       NaN    1.2002       Inf       Inf
 *
 * >> 2^(-24)
 * ans =
 * 5.9605e-008
 *
 * >> halfprecision(ans)
 * ans =
 *     1
 *
 * >> halfprecision(ans,'disp')
 * 5.9605e-008
 *
 * >> 2^(-25)
 * ans =
 * 2.9802e-008
 *
 * >> halfprecision(ans)
 * ans =
 *     1
 *
 * >> halfprecision(ans,'disp')
 * 5.9605e-008
 *
 * >> 2^(-26)
 * ans =
 *  1.4901e-008
 *
 * >> halfprecision(ans)
 * ans =
 *     0
 *
 * >> halfprecision(ans,'disp')
 *    0
 *
 * Note that the special cases of -Inf, +Inf, and NaN are handled correctly.
 * Also, note that the -1e30 and 1e30 values overflow the half precision format
 * and are converted into half precision -Inf and +Inf values, and stay that
 * way when they are converted back into doubles.
 *
 * For the denormalized cases, note that 2^(-24) is the smallest number that can
 * be represented in half precision exactly. 2^(-25) will convert to 2^(-24)
 * because of the rounding algorithm used, and 2^(-26) is too small and underflows
 * to zero.
 *
 ********************************************************************************/

// Includes -------------------------------------------------------------------


// #include "matlab_overload.h"
#include <math.h>
#include <stdio.h>

#include <stdint.h>
#include <string.h>
#include "mex.h"
#include <omp.h>

// Macros ---------------------------------------------------------------------

// Needed for older MATLAB versions that do not have mwSize
#ifdef DEFINEMWSIZE
#define mwSize int
#endif

#define INT16_TYPE int16_t
#define UINT16_TYPE uint16_t
#define INT32_TYPE int32_t
#define UINT32_TYPE uint32_t 

#define THREADS 28


// Global ---------------------------------------------------------------------

int next;  // Little Endian adjustment
int checkieee = 1;  // Flag to check for IEEE754, Endian, and word size

// Prototypes -----------------------------------------------------------------
template <typename dtype_out,typename dtype_in >
void singles2halfp(mxArray *target[], const mxArray *source[], const mwSize numel);


template <typename dtype_out,typename dtype_in >
void halfp2singles( mxArray *target[], const mxArray *source[], const mwSize numel); 


// overload the GetData function for each of the possible data type + select the correct matlab get function 
inline void  GetData(const mxArray *in, const mxComplexSingle *& out) {  out =  mxGetComplexSingles(in); return;}; 
inline void  GetData(const mxArray *in, const mxComplexUint16 *& out) {  out =  mxGetComplexUint16s(in); return;}; 
inline void  GetData(const mxArray *in, const mxSingle *& out)        {  out =  mxGetSingles(in); return;}; 
inline void  GetData(const mxArray *in, const mxUint16 *& out)        {  out =  mxGetUint16s(in); return;}; 
inline void  GetData(const mxArray *in, mxComplexSingle *& out)       {  out =  mxGetComplexSingles(in); return;}; 
inline void  GetData(const mxArray *in, mxComplexUint16 *& out)       {  out =  mxGetComplexUint16s(in); return;}; 
inline void  GetData(const mxArray *in, mxSingle *& out)              {  out =  mxGetSingles(in); return;}; 
inline void  GetData(const mxArray *in, mxUint16 *& out)              {  out =  mxGetUint16s(in); return;}; 


// Gateway Function -----------------------------------------------------------

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    mwSize ndim; // Number of dimensions of input
    mwSize numel; // Number of elements of input
    const mwSize *dims; // Pointer to dimensions array
    mxClassID classid; // Class id of input or desired output
    mxComplexity complexity; // Complexity of input
    mxArray *rhs[2], *lhs[1], *currentformat[1]; // Used for callbacks into MATLAB
    char *classname; // Class name of desired output
    int disp = 0; // Display flag
    double one = 1.0; // Used for checking IEEE754 floating point format
    UINT32_TYPE *ip; // Used for checking IEEE754 floating point format
    
    if( checkieee ) { // 1st call, so check for IEEE754, Endian, and word size
        ip = (UINT32_TYPE *) &one;
        if( *ip ) { // If Big Endian, then no adjustment
            next = 0;
        } else { // If Little Endian, then adjustment will be necessary
            next = 1;
            ip++;
        }
        if( *ip != 0x3FF00000u ) { // Check for exact IEEE 754 bit pattern of 1.0
            mexErrMsgTxt("Floating point bit pattern is not IEEE 754");
        }
        if( sizeof(INT16_TYPE) != 2 || sizeof(INT32_TYPE) != 4 ) {
            mexErrMsgTxt("Internal error. short is not 16-bits, or long is not 32-bits.");
        }
        checkieee = 0; // Everything checks out OK
    }

// Check input arguments for number and type
    
    if( nlhs > 1 ) {
        mexErrMsgTxt("Too many output arguments.");
    }
    if( nrhs != 1 && nrhs != 2 ) {
        mexErrMsgTxt("Need one or two input arguments.");
    }
    if( mxIsSparse(prhs[0]) ) {
        mexErrMsgTxt("Sparse matrices not supported.");
    }
    if( nrhs == 2 ) {
        if( !mxIsChar(prhs[1]) )
            mexErrMsgTxt("2nd input argument must be char string naming desired class, or 'disp'.");
        if( !mxIsInt16(prhs[0]) && !mxIsUint16(prhs[0]) )
            mexErrMsgTxt("1st input argument must be uint16 or int16 class.");
        classname = mxArrayToString(prhs[1]);
        if( strcmp(classname,"double") == 0 ) {  // Check 2nd input for proper class name string
            classid = mxDOUBLE_CLASS;
        } else if( strcmp(classname,"single") == 0 ) {
            classid = mxSINGLE_CLASS;
        } else if( strcmp(classname,"int8") == 0 ) {
            classid = mxINT8_CLASS;
        } else if( strcmp(classname,"uint8") == 0 ) {
            classid = mxUINT8_CLASS;
        } else if( strcmp(classname,"int16") == 0 ) {
            classid = mxINT16_CLASS;
        } else if( strcmp(classname,"uint16") == 0 ) {
            classid = mxUINT16_CLASS;
        } else if( strcmp(classname,"int32") == 0 ) {
            classid = mxINT32_CLASS;
        } else if( strcmp(classname,"uint32") == 0 ) {
            classid = mxUINT32_CLASS;
        } else if( strcmp(classname,"int64") == 0 ) {
            classid = mxINT64_CLASS;
        } else if( strcmp(classname,"uint64") == 0 ) {
            classid = mxUINT64_CLASS;
        } else if( strcmp(classname,"char") == 0 ) {
            classid = mxCHAR_CLASS;
        } else if( strcmp(classname,"logical") == 0 ) {
            classid = mxLOGICAL_CLASS;
        } else if( strcmp(classname,"disp") == 0 ) {
            disp = 1;
        } else {
            mexErrMsgTxt("2nd input argument must be char string naming desired class, or 'disp'.");
        }
    }
    
    if( nrhs == 1 ) {  // Convert from input MATLAB variable to halfprecision ---------------------
    
        classid = mxGetClassID(prhs[0]);  // Check for supported class
        switch( classid ) {
        case mxDOUBLE_CLASS:
        case mxSINGLE_CLASS:
        case mxUINT64_CLASS:
        case mxINT64_CLASS:
        case mxUINT32_CLASS:
        case mxINT32_CLASS:
        case mxUINT16_CLASS:
        case mxINT16_CLASS:
        case mxUINT8_CLASS:
        case mxINT8_CLASS:
        case mxLOGICAL_CLASS:
        case mxCHAR_CLASS:
            break;
        default:
            mexErrMsgTxt("Unable to convert input argument to halfprecision.");
        }
        complexity = mxIsComplex(prhs[0]) ? mxCOMPLEX : mxREAL; // Get stats of input variable
        numel = mxGetNumberOfElements(prhs[0]);
        ndim = mxGetNumberOfDimensions(prhs[0]);
        dims = mxGetDimensions(prhs[0]);
        plhs[0] = mxCreateNumericArray(ndim, dims, mxUINT16_CLASS, complexity); // halfprecision stored as uint16
        switch( classid ) {
        case mxDOUBLE_CLASS:  // Custom code for double class
            mexErrMsgTxt("Use single precision output");
            break;
        case mxSINGLE_CLASS:  // Custom code for single class
               if (mxIsComplex(prhs[0]))
                   singles2halfp<mxComplexUint16,mxComplexSingle>(plhs,prhs,2*numel);
               else
                   singles2halfp<mxUint16,mxSingle>(plhs,prhs,numel);
            break;
        case mxUINT64_CLASS:
        case mxINT64_CLASS:
        case mxUINT32_CLASS:
        case mxINT32_CLASS:
        case mxUINT16_CLASS:
        case mxINT16_CLASS:
        case mxUINT8_CLASS:
        case mxINT8_CLASS:
        case mxLOGICAL_CLASS:
        case mxCHAR_CLASS:  // All other classes get converted to single first
            mexErrMsgTxt("Use single precision output");
        }
    
    } else {  // Convert halfprecision to desired class ----------------------------------
        complexity = mxIsComplex(prhs[0]) ? mxCOMPLEX : mxREAL; // Get stats of input variable
        numel = mxGetNumberOfElements(prhs[0]);
        ndim = mxGetNumberOfDimensions(prhs[0]);
        dims = mxGetDimensions(prhs[0]);
        switch( classid ) {
        case mxDOUBLE_CLASS: // Custom code for double class
            mexErrMsgTxt("Use single precision input");
            break;
        case mxSINGLE_CLASS: // Custom code for single class
            plhs[0] = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, complexity);
            // Convert input to single
            if (mxIsComplex(prhs[0]))
               halfp2singles<mxComplexSingle, mxComplexUint16>(plhs,prhs,2*numel);
            else
               halfp2singles<mxSingle, mxUint16>(plhs,prhs,numel);

            break;
        case mxUINT64_CLASS:
        case mxINT64_CLASS:
        case mxUINT32_CLASS:
        case mxINT32_CLASS:
        case mxUINT16_CLASS:
        case mxINT16_CLASS:
        case mxUINT8_CLASS:
        case mxINT8_CLASS:
        case mxLOGICAL_CLASS:
        case mxCHAR_CLASS: // All other classes get converted to single first
            mexErrMsgTxt("Use single precision input");
            break;
        default:
            mexErrMsgTxt("Unable to convert input argument to halfprecision.");
        }
        mxFree(classname);
    }
    
}

//-----------------------------------------------------------------------------

template <typename dtype_out,typename dtype_in >
void singles2halfp( mxArray *target[], const mxArray *source[], const mwSize numel)
{
    
    
    dtype_in const * source_array;
    GetData(source[0], source_array);
    dtype_out * target_array;
    GetData(target[0], target_array);
 

    UINT16_TYPE *hp = (UINT16_TYPE *) target_array; // Type pun output as an unsigned 16-bit int
    UINT32_TYPE const *xp = (UINT32_TYPE *) source_array; // Type pun input as an unsigned 32-bit int
   
    UINT16_TYPE    hs, he, hm;
    UINT32_TYPE x, xs, xe, xm;
    int hes;
    mwSize i;
    
    if( source == NULL || target == NULL ) { // Nothing to convert (e.g., imag part of pure real)
        return;
    }
    
    #pragma omp parallel for schedule(static) num_threads(THREADS) private(i,x,xs,xe,xm,hs,he,hm,hes)
    for( i = 0; i < numel; i++)
    {
        x = xp[i];
        if( (x & 0x7FFFFFFFu) == 0 ) {  // Signed zero
            hp[i] = (UINT16_TYPE) (x >> 16);  // Return the signed zero
        } else { // Not zero
            xs = x & 0x80000000u;  // Pick off sign bit
            xe = x & 0x7F800000u;  // Pick off exponent bits
            xm = x & 0x007FFFFFu;  // Pick off mantissa bits
            if( xe == 0 ) {  // Denormal will underflow, return a signed zero
                hp[i] = (UINT16_TYPE) (xs >> 16);
            } else if( xe == 0x7F800000u ) {  // Inf or NaN (all the exponent bits are set)
                if( xm == 0 ) { // If mantissa is zero ...
                    hp[i] = (UINT16_TYPE) ((xs >> 16) | 0x7C00u); // Signed Inf
                } else {
                    hp[i] = (UINT16_TYPE) 0xFE00u; // NaN, only 1st mantissa bit set
                }
            } else { // Normalized number
                hs = (UINT16_TYPE) (xs >> 16); // Sign bit
                hes = ((int)(xe >> 23)) - 127 + 15; // Exponent unbias the single, then bias the halfp
                if( hes >= 0x1F ) {  // Overflow
                    hp[i] = (UINT16_TYPE) ((xs >> 16) | 0x7C00u); // Signed Inf
                } else if( hes <= 0 ) {  // Underflow
                    if( (14 - hes) > 24 ) {  // Mantissa shifted all the way off & no rounding possibility
                        hm = (UINT16_TYPE) 0u;  // Set mantissa to zero
                    } else {
                        xm |= 0x00800000u;  // Add the hidden leading bit
                        hm = (UINT16_TYPE) (xm >> (14 - hes)); // Mantissa
                        if( (xm >> (13 - hes)) & 0x00000001u ) // Check for rounding
                            hm += (UINT16_TYPE) 1u; // Round, might overflow into exp bit, but this is OK
                    }
                    hp[i] = (hs | hm); // Combine sign bit and mantissa bits, biased exponent is zero
                } else {
                    he = (UINT16_TYPE) (hes << 10); // Exponent
                    hm = (UINT16_TYPE) (xm >> 13); // Mantissa
                    if( xm & 0x00001000u ) // Check for rounding
                        hp[i] = (hs | he | hm) + (UINT16_TYPE) 1u; // Round, might overflow to inf, this is OK
                    else
                        hp[i] = (hs | he | hm);  // No rounding
                }
            }
        }
    }  
}

template <typename dtype_out,typename dtype_in>
void halfp2singles( mxArray *target[], const mxArray *source[], const mwSize numel)
{
    
    dtype_in  * source_array;
    GetData(source[0], source_array);
    dtype_out const * target_array;
    GetData(target[0], target_array);
 

    UINT16_TYPE const *hp = (UINT16_TYPE *) source_array; // Type pun output as an unsigned 16-bit int
    UINT32_TYPE  *xp = (UINT32_TYPE *) target_array; // Type pun input as an unsigned 32-bit int
     
    UINT16_TYPE h, hs, he, hm;
    UINT32_TYPE xs, xe, xm;
    INT32_TYPE xes;
    int e;
    mwSize i; 
    
    if( source == NULL || target == NULL ) // Nothing to convert (e.g., imag part of pure real)
        return;
    
    #pragma omp parallel for schedule(static) num_threads(THREADS) private(i,e,xs,xe,xm,h,hs,he,hm,xes)
    for( i = 0; i < numel; i++)
    {
        h = hp[i];
        if( (h & 0x7FFFu) == 0 ) {  // Signed zero
            xp[i] = ((UINT32_TYPE) h) << 16;  // Return the signed zero
        } else { // Not zero
            hs = h & 0x8000u;  // Pick off sign bit
            he = h & 0x7C00u;  // Pick off exponent bits
            hm = h & 0x03FFu;  // Pick off mantissa bits
            if( he == 0 ) {  // Denormal will convert to normalized
                e = -1; // The following loop figures out how much extra to adjust the exponent
                do {
                    e++;
                    hm <<= 1;
                } while( (hm & 0x0400u) == 0 ); // Shift until leading bit overflows into exponent bit
                xs = ((UINT32_TYPE) hs) << 16; // Sign bit
                xes = ((INT32_TYPE) (he >> 10)) - 15 + 127 - e; // Exponent unbias the halfp, then bias the single
                xe = (UINT32_TYPE) (xes << 23); // Exponent
                xm = ((UINT32_TYPE) (hm & 0x03FFu)) << 13; // Mantissa
                xp[i] = (xs | xe | xm); // Combine sign bit, exponent bits, and mantissa bits
            } else if( he == 0x7C00u ) {  // Inf or NaN (all the exponent bits are set)
                if( hm == 0 ) { // If mantissa is zero ...
                    xp[i] = (((UINT32_TYPE) hs) << 16) | ((UINT32_TYPE) 0x7F800000u); // Signed Inf
                } else {
                    xp[i] = (UINT32_TYPE) 0xFFC00000u; // NaN, only 1st mantissa bit set
                }
            } else { // Normalized number
                xs = ((UINT32_TYPE) hs) << 16; // Sign bit
                xes = ((INT32_TYPE) (he >> 10)) - 15 + 127; // Exponent unbias the halfp, then bias the single
                xe = (UINT32_TYPE) (xes << 23); // Exponent
                xm = ((UINT32_TYPE) hm) << 13; // Mantissa
                xp[i] = (xs | xe | xm); // Combine sign bit, exponent bits, and mantissa bits
            }
        }
    }
}

