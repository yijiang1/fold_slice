/*
 * Copyright (c) 2010,2011 Joshua V Dillon
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or
 * without modification, are permitted provided that the
 * following conditions are met:
 *  * Redistributions of source code must retain the above
 *    copyright notice, this list of conditions and the
 *    following disclaimer.
 *  * Redistributions in binary form must reproduce the above
 *    copyright notice, this list of conditions and the
 *    following disclaimer in the documentation and/or other
 *    materials provided with the distribution.
 *  * Neither the name of the author nor the names of its
 *    contributors may be used to endorse or promote products
 *    derived from this software without specific prior written
 *    permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY JOSHUA V DILLON ''AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL JOSHUA
 * V DILLON BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef sharedmatrix_hxx
#define sharedmatrix_hxx

/* Uses mxArray_tag as defined in matrix.h */
/*#define COMPLIANTMODE*/

/* Verbose outputs */
// #define DEBUG

/* This copies memory across, defeating the purpose of this function but useful for testing */ 
// #define SAFEMODE

/*
 * sharedmatrix.h
 *
 * This Matlab Mex program allows you to share matlab cell arrays and matrices
 * between different Maatlab processes.  It has four functions:
 * 1) Clone:
 *    Serialize a 2D cell array of 2D non-/sparse matrices or a 2D
 *    non-/sparse matrix to shared memory.
 * 2) Attach:
 *    "Reconstitute" the shared data into the appropriate Matlab object using
 *    shallow copying.
 * 3) Detach:
 *    Remove the shallow references and detach the shared memory from the
 *    Matlab data space.
 * 4) Free:
 *    Mark the shared memory for destruction.
 *
 * Written by Joshua V Dillon
 * August 27, 2010
 *
 * Additional Contributors"
 * Andrew Smith - Februrary 18, 2011
 *
 * Revision History:
 * Sep. 9. 2010  Corrected the mishandling of sparse logical matrices.
 * Apr. 8, 2011  - Merged Andrew Smith's struct code
 *               - Added Andrew Smith's windows/boost contribution
 */

/*
 * This code can be compiled from within Matlab or command-line, assuming the
 * system is appropriately setup.  To compile, invoke:
 *
 * For 32-bit machines:
 *     mex -O -v sharedmatrix.c
 * For 64-bit machines:
 *     mex -largeArrayDims -O -v sharedmatrix.c
 *
 */

/* 
 * Programmer's Notes:
 *
 * MEX C API:
 * http://www.mathworks.com/access/helpdesk/help/techdoc/apiref/bqoqnz0.html
 *
 * Testing:
 *
x=sparse((rand(3,4)<.5).*rand(3,4));
x=cell(2,2);x{1}=rand(3,4);x{4}=sparse((rand(2,4)>.5).*rand(2,4));
shmsiz=sharedmatrix('clone',12345,x)
y=sharedmatrix('attach',12345)
sharedmatrix('detach',12345,y)
 *
 */

/*
 * Possibly useful information/related work:
 *
 * http://www.mathworks.com/matlabcentral/fileexchange/24576
 * http://www.mathworks.in/matlabcentral/newsreader/view_thread/254813
 * http://www.mathworks.com/matlabcentral/newsreader/view_thread/247881#639280
 * 
 * http://groups.google.com/group/comp.soft-sys.matlab/browse_thread/thread/c241d8821fb90275/47189d498d1f45b8?lnk=st&q=&rnum=1&hl=en#47189d498d1f45b8
 * http://www.mk.tu-berlin.de/Members/Benjamin/mex_sharedArrays
 */

/* Possibily useful undocumented functions (see links at end for details): */
/* extern mxArray *mxCreateSharedDataCopy(const mxArray *pr);              */
/* extern bool mxUnshareArray(const mxArray *pr, const bool noDeepCopy);   */
/* extern mxArray *mxUnreference(const mxArray *pr);                       */

#ifdef DEBUG
	#ifdef __STDC__
		/*#define DEBUG_HEADER "[%s:% 4d] % 15s(): "*/
		#define DEBUG_HEADER  "\033[1;34;40m[%s:% 4d] \033[1;31;40m% 15s(): \033[0m"
		#define mydebug(format,args...)               \
			mexPrintf(DEBUG_HEADER format "\n",       \
				__FILE__,__LINE__,__FUNCTION__,##args)
	#else
		#define DEBUG_HEADER "DEBUG: "
		#define mydebug(format,args...) \
			mexPrintf(DEBUG_HEADER format "\n",##args)
	#endif
#else
	#define mydebug(format,args...) (void)0
#endif

#ifdef COMPLIANTMODE

/* We will be accessing the mxArray parts directly via the mxArray_tag struct
 * defined in matrix.h.  */
#define ARRAY_ACCESS_INLINING
#include "matrix.h"

#else

struct mxArray_tag {
	void    *reserved;
	int      reserved1[2];
	void    *reserved2;
	size_t  number_of_dims;
	unsigned int reserved3;
	struct {
		unsigned int    flag0 : 1;
		unsigned int    flag1 : 1;
		unsigned int    flag2 : 1;
		unsigned int    flag3 : 1;
		unsigned int    flag4 : 1;
		unsigned int    flag5 : 1;
		unsigned int    flag6 : 1;
		unsigned int    flag7 : 1;
		unsigned int    flag7a: 1;
		unsigned int    flag8 : 1;
		unsigned int    flag9 : 1;
		unsigned int    flag10 : 1;
		unsigned int    flag11 : 4;
		unsigned int    flag12 : 8;
		unsigned int    flag13 : 8;
	}   flags;
	size_t reserved4[2];
	union {
		struct {
			void  *pdata;
			void  *pimag_data;
			void  *reserved5;
			size_t reserved6[3];
		}   number_array;
	}   data;
};

#endif
typedef struct mxArray_tag mxArrayHack;

/* standard mex include; after hack */
#include "mex.h"

/* max length of directive string */
#define MAXDIRECTIVELEN 256

/* these are used for recording structure field names */
const char term_char = ';';     /*use this character to terminate a string containing the list of fields.  Do this because it can't be in a valid field name*/
const size_t  align_size = 8;   /*the pointer alignment size, so if pdata is a valid pointer then &pdata[i*align_size] will also be.  Ensure this is >= 4*/

/* 
 * The header_t object will be copied to shared memory in its entirety.
 *
 * Immediately after each copied header_t will be the matrix data values 
 * [size array, field names, (_at_most_ the four arrays [pr,pi,ir,jc] and in this order)].
 *
 * The data_t objects will never be copied to shared memory and serve only 
 * to abstract away mex calls and simplify the deep traversals in matlab. 
 *
 */

typedef struct data data_t;
typedef struct header header_t;

/* structure used to record all of the data addresses */
struct data {
	mwSize    *pSize;      /* pointer to the size array */   
	void*      pr;         /* real data portion */
	void*      pi;         /* imaginary data portion */
	mwIndex   *ir;         /* row indexes, for sparse */
	                       /* OR: may also be list of a structures fields,
	                              each field name will be seperated by a null
	                              character and terminated with a ";" */
	mwIndex   *jc;         /* cumulative column counts, for sparse */
	data_t    *child_dat;  /* array of children data structures, for cell */
	header_t  *child_hdr;  /* array of corresponding children header structures, for cell */
};

/* captures fundamentals of the mxArray */
/* In the shared memory the storage order is [header, size array, field_names,
 * real dat, image data, sparse index r, sparse index c]  */
struct header {
	bool       isCell;
	bool       isSparse;
	bool       isComplex;
	bool       isStruct;      
	mxClassID  classid;  /* matlab class id */
	size_t     nDims;    /* dimensionality of the matrix.  The size array immediately follows the header */ 
	size_t     elemsiz;  /* size of each element in pr and pi */
	size_t     nzmax;    /* length of pr,pi */
	size_t     nFields;  /* the number of fields.  The field string immediately follows the size array */
	size_t     shmsiz;   /* size of serialized object (header + size array + field names string) */
};

/* Remove shared memory references to input matrix (in-situ), recursively    */
/* if needed.                                                                */
char*  deepdetach     (mxArray *mxInput);

/* Shallow copy matrix from shared memory into Matlab form.                  */
size_t shallowrestore (char *shm, mxArray** p_mxInput);

/* Recursively descend through Matlab matrix to assess how much space its    */
/* serialization will require.                                               */
size_t deepscan       (header_t *hdr, data_t *dat, const mxArray* mxInput);

/* Descend through header and data structure and copy relevent data to       */
/* shared memory.                                                            */
void   deepcopy       (header_t *hdr, data_t *dat, char *shared_mem, bool allocate_only);

/* Descend through header and data structure and free the memory.            */
void   deepfree       (data_t *dat);

/* Pads the size to something that guarantees pointer alignment.             */
__inline size_t pad_to_align(size_t size) {
	if (size % align_size)
		size += align_size - (size % align_size);
	return size;
}

/* Function to find the number of bytes required to store all of the         */
/* field names of a structure                                                */
int FieldNamesSize(const mxArray * mxStruct);

/* Function to copy all of the field names to a character array              */
/* Use FieldNamesSize() to allocate the required size of the array           */
/* returns the number of bytes used in pList                                 */
int CopyFieldNames(const mxArray * mxStruct, char* pList);

/* This function finds the number of fields contained within in a string     */
/* the string is terminated by term_char, a character that can't be in a     */
/* field name.  pBytes is always an aligned number                           */
int NumFieldsFromString(const char* pString, size_t *pfields, size_t* pBytes);

/* Function to take point a each element of the char** array at a list of    */
/* names contained in string                                                 */
/* ppCharArray must be allocated to length num_names                         */
/* names are seperated by null termination characters                        */
/* each name must start on an aligned address (see CopyFieldNames())         */
/* e.g. pCharArray[0] = name_1, pCharArray[1] = name_2 ...                   */
/* returns 0 if successful                                                   */
int PointCharArrayAtString(char ** pCharArray, char* pString, int nFields);

/* Function to find the bytes in the string starting from the end of the     */
/* string returns < 0 on error                                               */
int BytesFromStringEnd(const char* pString, size_t* pBytes);

// #ifdef SAFEMODE
/* A convenient function for safe assignment of memory to an mxArray         */
void* safeCopy(void* pBuffer, mwSize Bytes) {
	void* pSafeBuffer;
    
	/* ensure Matlab knows it */
	pSafeBuffer = mxMalloc(Bytes);
	if (pSafeBuffer != NULL)
		memcpy(pSafeBuffer, pBuffer, Bytes); /* and copy the data */

	return pSafeBuffer;
}
// #endif

#endif

