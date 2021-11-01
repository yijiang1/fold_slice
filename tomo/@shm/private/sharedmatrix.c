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


#include <string.h>
#include <sys/shm.h>
#include <ctype.h>
#include "sharedmatrix.h"


/* ------------------------------------------------------------------------- */
/* Matlab gateway function                                                   */
/*                                                                           */
/* (see sharedmatrix.m for description)                                      */
/* ------------------------------------------------------------------------- */
void mexFunction( int nlhs,       mxArray *plhs[], 
                  int nrhs, const mxArray *prhs[]  )
{
	/* mex inputs */
	const mxArray  *mxDirective; /*  directive {clone, attach, detach, free} */
	const mxArray  *mxShmKey;    /*  name of the shared memory segment       */ 
	const mxArray  *mxInput;     /*  input array (for clone)                 */
	/* mex outputs */
	mxArray        *mxOutput=NULL;

	/* for storing directive (string) input */
	char directive[MAXDIRECTIVELEN+1];

	/* for working with shared memory ... */
	key_t     shmkey=0;
	size_t    shmsiz=0;     /* apparently you can attach w/o specifying this */
	size_t    shmsiz_tmp=0; /* for shallowrestore recursion */
	int       shmid;
	char     *shm=NULL;

	/* for storing the mxArrays ... */
	header_t  hdr;
	data_t    dat;

	/* check min number of arguments */
	if ( nrhs<2 )
		mexErrMsgIdAndTxt("MATLAB:sharedmatrix",
			"Minimum input arguments missing; must supply directive and key.");

	/* assign inputs */
	mxDirective = prhs[0];
	mxShmKey = prhs[1];

	/* get directive (ARGUMENT 0) */
	if ( mxGetString(mxDirective,directive,MAXDIRECTIVELEN)!=0 )
		mexErrMsgIdAndTxt("MATLAB:sharedmatrix",
			"First input argument must be {'clone','attach','detach','free'}.");

	/* get key (ARGUMENT 1) */
	if ( mxIsNumeric(mxShmKey) && mxGetNumberOfElements(mxShmKey)==1 )
		shmkey = (key_t)mxGetScalar(mxShmKey);
	else
		mexErrMsgIdAndTxt("MATLAB:sharedmatrix",
			"Second input argument must be a valid key (numeric scalar).");

	/* check outputs */
	if ( nlhs > 1 )
		mexErrMsgIdAndTxt("MATLAB:sharedmatrix",
			"Function returns only one value.");

	/* clone, attach, detach, free */
	switch ( tolower(directive[0]) ) {

	/* --------------------------------------------------------------------- */
	/*    Clone                                                              */
	/* --------------------------------------------------------------------- */
	case 'c':

		if ( nrhs<3 )
			mexErrMsgIdAndTxt("MATLAB:sharedmatrix:clone",
				"Required third argument missing (variable).");

		if ( mxGetNumberOfElements(prhs[2])<2 )
			mexErrMsgIdAndTxt("MATLAB:sharedmatrix:clone",
				"Required third argument (variable) must have at least two elements.");

		/* Assign (ARGUMENT 2) */
		mxInput = prhs[2];
        
        bool allocate_only = false; 
        if (nrhs>3) {
            int * tmp = (int*)mxGetData(prhs[3]);
            allocate_only = tmp[0];
        }

		/* scan input data */
		shmsiz = deepscan(&hdr,&dat,mxInput);
		mydebug("clone: deepscan done (%u)",(unsigned)shmsiz);
		
		/* create the segment */
		if ( (shmid=shmget(shmkey,shmsiz,IPC_CREAT|IPC_EXCL|0644 )) < 0 )
			mexErrMsgIdAndTxt("MATLAB:sharedmatrix:clone",
				"Unable to create shared memory segment.");
		
		/* attach the segment to this dataspace */
		shm = (char*)shmat(shmid,NULL,0);
		if ( shm==(char*)-1 )
			mexErrMsgIdAndTxt("MATLAB:sharedmatrix:clone",
				"Unable to attach shared memory to data space.");

		mydebug("clone: shared memory at: 0x%llx",shm);

		/* copy data to the shared memory */
		deepcopy(&hdr,&dat,shm, allocate_only);
		mydebug("clone: deep copy successful");

		/* free temporary allocation */
		deepfree(&dat);
		shmdt(shm);

		/* return its size */
		mxOutput = mxCreateDoubleScalar((double)shmsiz);

		break;

	/* --------------------------------------------------------------------- */
	/*    Attach                                                             */
	/* --------------------------------------------------------------------- */
	case 'a':

		/* get SM id */
		if ( (shmid=shmget(shmkey,shmsiz,IPC_CREAT|0644)) < 0 )
			mexErrMsgIdAndTxt("MATLAB:sharedmatrix:attach",
				"Unable to get shared memory id.");

		/* attach the segment to this dataspace */
		shm = (char*)shmat(shmid,NULL,0);
		if ( shm==(char*)-1 )
			mexErrMsgIdAndTxt("MATLAB:sharedmatrix:attach",
				"Unable to attach shared memory to data space.");

		mydebug("attach: shared memory at: 0x%llx",shm);
		
		/* restore the segments (ignore return value) */
		shmsiz_tmp=shallowrestore(shm,&mxOutput);

		mydebug("attach: done");

		break;

	/* --------------------------------------------------------------------- */
	/*    Detach                                                             */
	/* --------------------------------------------------------------------- */
	case 'd':

		if ( nrhs<3 || mxIsEmpty(prhs[2]) )
			mexErrMsgIdAndTxt("MATLAB:sharedmatrix:detach",
				"Required third argument (variable) missing or empty.");

		/* Assign */
		mxInput = prhs[2];

		shm = deepdetach((mxArray*)mxInput); /* Abuse const */
		mxOutput = mxCreateDoubleMatrix(0,0,mxREAL);

		mydebug("detach: shared memory at: 0x%llx",shm);

		if ( shmdt(shm) != 0 )
			mexErrMsgIdAndTxt("MATLAB:sharedmatrix:detach",
				"Unable to detach shared memory.");

		break;

	/* --------------------------------------------------------------------- */
	/*    Free                                                               */
	/* --------------------------------------------------------------------- */
	case 'f':

		/* get SM id */
		if ( (shmid=shmget(shmkey,shmsiz,0644)) < 0 )
			mexErrMsgIdAndTxt("MATLAB:sharedmatrix:free",
				"Unable to get shared memory id.");

		/* free SM */
		if ( (shmctl(shmid,IPC_RMID,(struct shmid_ds *)NULL)) != 0 )
			mexErrMsgIdAndTxt("MATLAB:sharedmatrix:free",
				"Unable to destroy shared memory.");

		mxOutput = mxCreateDoubleMatrix(0,0,mxREAL);

		break;

	/* --------------------------------------------------------------------- */
	/*    Error                                                              */
	/* --------------------------------------------------------------------- */
	default:
		mexErrMsgIdAndTxt("MATLAB:sharedmatrix",
			"Unrecognized directive.");
	}

	/* assign the output */
	plhs[0] = mxOutput;

}




/* ------------------------------------------------------------------------- */
/* deepdetach                                                                */
/*                                                                           */
/* Remove shared memory references to input matrix (in-situ), recursively    */
/* if needed.                                                                */
/*                                                                           */
/* Arguments:                                                                */
/*    Input matrix to remove references.                                     */
/* Returns:                                                                  */
/*    Pointer to start of shared memory segment.                             */
/* ------------------------------------------------------------------------- */
char* deepdetach(mxArray *mxInput) {

	/* uses side-effects! */
	mwSize dims[] = {0,0};
	mwSize nzmax  = 0;
	size_t elemsiz;
	size_t i,j,n,m;
	size_t offset=0;
	char *shm=(char*)NULL,*tmp=(char*)NULL;


	/* restore matlab memory */
	if ( mxInput==(mxArray*)NULL || mxIsEmpty(mxInput) ) {
		/* do nothing */

	} else if ( mxIsStruct(mxInput) ) {
		/* struct case */

		/* detach each field for each element */
		m = mxGetNumberOfElements(mxInput);
		n = mxGetNumberOfFields(mxInput);
		for ( i=0; i<m; i++ ) {			/* element */
			for ( j=0; j<n; j++ ) {		/* field */
				/* detach this one */
				tmp=deepdetach((mxArray*)(mxGetFieldByNumber(mxInput,i,j)));
				if ( shm==NULL ) {
					if ( tmp!=NULL ) shm=tmp;
					offset += FieldNamesSize(mxInput); /* always a multiple of alignment size */
					offset += sizeof(header_t);
					offset += pad_to_align( mxGetNumberOfDimensions(mxInput) * sizeof(mwSize) );
					/*offset += strBytes;*/
				}
			}
		}
		shm -= offset;

	} else if ( mxIsCell(mxInput) ) {
		/* cell case */

		/* detach each element */
		m = mxGetNumberOfElements(mxInput);
		for ( i=0; i<m; i++ ) {
			/* detach this one */
			tmp=deepdetach((mxArray*)(mxGetCell(mxInput,i)));
			if ( shm==NULL ) {
				if ( tmp!=NULL ) shm=tmp;
				offset += sizeof(header_t);
				offset += pad_to_align( mxGetNumberOfDimensions(mxInput) * sizeof(mwSize) );
			}
		}
		shm -= offset;

	} else if ( mxIsNumeric(mxInput) || mxIsLogical(mxInput) || mxIsChar(mxInput)) {
		/* matrix case */
		
		n = sizeof(header_t);
		n += pad_to_align( mxGetNumberOfDimensions(mxInput) * sizeof(mwSize) );
		shm = (char*)(mxGetData(mxInput) - n);

		/* in safe mode these entries were allocated so remove them properly */
		if (  mxIsComplex(mxInput) )  {
            mxFree(mxGetData(mxInput));
            mxFree(mxGetImagData(mxInput));
        }

		/* handle sparse objects */
// 		if ( mxIsSparse(mxInput) ) {
// 			/* can't give sparse arrays zero size (nzmax must be 1) */
// 			dims[0] = dims[1] = nzmax = 1;
// 			if ( mxSetDimensions(mxInput,dims,2) )
// 				mexErrMsgIdAndTxt("MATLAB:sharedmatrix:detach",
// 					"Unable to resize the array.");
// 
// 			/* in safe mode these entries were allocated so remove them
// 			 * properly */
// 			#ifdef SAFEMODE
// 			mxFree(mxGetIr(mxInput));
// 			mxFree(mxGetJc(mxInput));
// 			#endif
// 
// 			/* allocate 1 element */
// 			elemsiz = mxGetElementSize(mxInput);
// 			mxSetData(mxInput,mxCalloc(nzmax,elemsiz));
// 			if ( mxIsComplex(mxInput) )
// 				mxSetImagData(mxInput,mxCalloc(nzmax,elemsiz));
// 			else
// 				mxSetImagData(mxInput,(void*)NULL);
// 			
// 			/* allocate 1 element */
// 			mxSetNzmax(mxInput,nzmax);
// 			mxSetIr(mxInput,(mwSize*)mxCalloc(nzmax,sizeof(mwIndex)));
// 			mxSetJc(mxInput,(mwSize*)mxCalloc(dims[1]+1,sizeof(mwIndex)));
// 
// 		} else {
			/* Can have zero size, so nullify data storage containers */
			if (mxSetDimensions(mxInput,dims,2))
				mexErrMsgIdAndTxt("MATLAB:sharedmatrix:detach",
					"Unable to resize the array.");

			/* Doesn't allocate or deallocate any space for the pr or pi arrays */
			mxSetData(mxInput, (void*)NULL);
			mxSetImagData(mxInput, (void*)NULL);
// 		}

	} else {
		mexErrMsgIdAndTxt("MATLAB:sharedmatrix:detach",
			"Unsupported type.");

	}

	return shm;
}



/* ------------------------------------------------------------------------- */
/* shallowrestore                                                            */
/*                                                                           */
/* Shallow copy matrix from shared memory into Matlab form.                  */
/*                                                                           */
/* Arguments:                                                                */
/*    shared memory segment.                                                 */
/*    pointer to an mxArray pointer.                                         */
/* Returns:                                                                  */
/*    size of shared memory segment.                                         */
/* ------------------------------------------------------------------------- */
size_t shallowrestore(char *shm, mxArray** p_mxInput)  {	

	/* for working with shared memory ... */
	size_t i,shmsiz;
	mxArray *mxChild;	/* temporary pointer */

	/* for working with payload ... */
	header_t *hdr;
	mwSize   *pSize;             /* pointer to the array dimensions */
	void     *pr=NULL,*pi=NULL;  /* real and imaginary data pointers */
	mwIndex  *ir=NULL,*jc=NULL;  /* sparse matrix data indices */

	/* for structures */
	int ret;				     
	size_t nFields;		/* number of fields */
	size_t field_num;	/* current field */
	size_t strBytes;	/* Number of bytes in the string recording the field names */
	char **pFieldStr;	/* Array to pass to Matlab with the field names pFields[0] is a poitner to the first field name, pFields[1] is the pointer to the second*/

	/* retrieve the data */
	hdr = (header_t*)shm;
	shm += sizeof(header_t);

	/* the size pointer */
	pSize = (mwSize*)shm;

	/* skip over the stored sizes */
	shm += pad_to_align(sizeof(mwSize)*hdr->nDims);

	if ( hdr->isStruct ) {
		/* struct case */

		mydebug("shallowrestore: found structure %d x %d",pSize[0],pSize[1]);

		/* Pull out the field names, they follow the size*/
		ret = NumFieldsFromString(shm, &nFields, &strBytes);
		
		/* check the recovery */
		if ((ret) || (nFields != hdr->nFields))
			mexErrMsgIdAndTxt("MATLAB:sharedmatrix:attach",
				"Structure fields have not been recovered properly.");

		/* And convert to something matlab understands */
		pFieldStr = (char**)mxCalloc(nFields, sizeof(char*));
		PointCharArrayAtString(pFieldStr, shm, nFields);
		
		/* skip over the stored field string */
		shm += strBytes;

		/*create the matrix */
		*p_mxInput = mxCreateStructArray(hdr->nDims, pSize, hdr->nFields, (const char**)pFieldStr);
		mxFree(pFieldStr);  /* no longer need it */

		for ( i=0; i<hdr->nzmax; i++ ) {	/* each element */
			for ( field_num=0; field_num<hdr->nFields; field_num++ ) { /* each field */
				mydebug("shallowrestore: working on %d field %d at 0x%llx",i,field_num, *p_mxInput);

				/* And fill it */
				shmsiz = shallowrestore(shm,&mxChild);					  /* restore the mxArray */
				mxSetFieldByNumber(*p_mxInput, i, field_num, mxChild);    /* and pop it in	   */
				shm += shmsiz;
			
				mydebug("shallowrestore: completed %d field %d",i, field_num);
			}
		}

	} else if ( hdr->isCell ) {
		/* cell case */

		mydebug("shallowrestore: found cell %d x %d",pSize[0], pSize[1]);

		/* Create the array */
		*p_mxInput = mxCreateCellArray(hdr->nDims, pSize);
		for ( i=0; i<hdr->nzmax; i++ ) {
			mydebug("shallowrestore: working on %d at 0x%llx",i,mxChild);
			
			/* And fill it */
			shmsiz = shallowrestore(shm,&mxChild);
			mxSetCell(*p_mxInput, i, mxChild);
			shm += shmsiz;
			
			mydebug("shallowrestore: completed %d at 0x%llx",i,mxChild);
		}


	} else {	
		/* matrix case */

		/* this is the address of the first data */ 
		pr = (void*)shm;
		shm += pad_to_align((hdr->nzmax)*(hdr->elemsiz));		/* takes us to the end of the real data */

 		/* if complex get a pointer to the complex data */
 		if ( hdr->isComplex ) {
 			pi = (void*)shm;
 			shm += pad_to_align((hdr->nzmax)*(hdr->elemsiz));	/* takes us to the end of the complex data */
 		}
		
		/* if sparse get a list of the elements */
// 		if ( hdr->isSparse ) {
// 			mydebug("shallowrestore: found non-cell, sparse 0x%llx",*p_mxInput);
// 
// 			ir = (mwIndex*)shm;
// 			shm += pad_to_align((hdr->nzmax)*sizeof(mwIndex));
// 
// 			jc = (mwIndex*)shm;
// 			shm += pad_to_align((pSize[1]+1)*sizeof(mwIndex));
// 
// 			/* need to create sparse logical differently */
// 			if ( hdr->classid==mxLOGICAL_CLASS )
// 				*p_mxInput = mxCreateSparseLogicalMatrix(0,0,0);
// 			else
// 				*p_mxInput = mxCreateSparse(0,0,0,(hdr->isComplex)?mxCOMPLEX:mxREAL);
// 
// 			/* free the memory it was created with */
// 			mxFree(mxGetIr(*p_mxInput));
// 			mxFree(mxGetJc(*p_mxInput));
// 
// 			/* set the size */
// 			if (mxSetDimensions(*p_mxInput, pSize, hdr->nDims)) 
// 				mexErrMsgIdAndTxt("MATLAB:sharedmatrix:attach",
// 					"Unable to resize the sparse array.");
// 
// 			/* set the pointers relating to sparse (set the real and imaginary data later)*/
// 			#ifndef SAFEMODE
// 		
// 			/* Hack  Method */
// 			((mxArrayHack*)*p_mxInput)->data.number_array.reserved5    = ir;
// 			((mxArrayHack*)*p_mxInput)->data.number_array.reserved6[0] = (size_t)jc;
// 			((mxArrayHack*)*p_mxInput)->data.number_array.reserved6[1] = (size_t)(hdr->nzmax);
// 
// 			/* dont do this as it makes a write turn into a copy to local memory */
// 			/*((mxArrayHack*)*p_mxInput)->reserved3 = 1;*/ /*give it a reference count*/
// 
// 			#else
// 		
// 			/* Safe - but takes it out of global memory */
// 			mxSetIr(*p_mxInput, (mwIndex*) safeCopy(ir,(hdr->nzmax)*sizeof(mwIndex)));
// 			mxSetJc(*p_mxInput, (mwIndex*) safeCopy(jc,(pSize[1]+1)*sizeof(mwIndex)));  						
// 			
// 			#endif
// 			mxSetNzmax(*p_mxInput,hdr->nzmax);
// 
// 		} else {

			/*  rebuild constitute the array  */
			mydebug("shallowrestore: found non-cell, non-sparse");

			/*  create an empty array */
			*p_mxInput = mxCreateNumericMatrix(0,0,hdr->classid,(hdr->isComplex)?mxCOMPLEX:mxREAL);

			mydebug("%d %d (of total: %d)",pSize[0],pSize[1],hdr->nDims);

			/* set the size */
			if (mxSetDimensions(*p_mxInput, pSize, hdr->nDims))
				mexErrMsgIdAndTxt("MATLAB:sharedmatrix:attach",
					"Unable to resize the array.");
// 		}

		
		/* Free the memory it was created with (shouldn't be necessary) */
		mxFree(mxGetPr(*p_mxInput));
		mxFree(mxGetPi(*p_mxInput));

		/* set the real and imaginary data */
		if ( ! hdr->isComplex ) 
        {
		/* Hack  Method */
		
		mydebug("pr:0x%llx",pr);

		/* BUGFIX: this line causes matlab ABORT when shared object is a scalar! */
		/* the theory is that scalars are copied not copy on reference hence
		 * the return of the function causes a failure. */
		((mxArrayHack*)(*p_mxInput))->data.number_array.pdata = pr;

 		if ( hdr->isComplex ) 
 			((mxArrayHack*)(*p_mxInput))->data.number_array.pimag_data = pi;
        }
		else
        {
		/* Safe - but takes it out of global memory */
        

         mxSetData(*p_mxInput, safeCopy(pr, hdr->nzmax*hdr->elemsiz));
		if ( hdr->isComplex )
			mxSetImagData(*p_mxInput, safeCopy(pi, hdr->nzmax*hdr->elemsiz));
        }

	}
	return hdr->shmsiz;
}



/* ------------------------------------------------------------------------- */
/* deepscan                                                                  */
/*                                                                           */
/* Recursively descend through Matlab matrix to assess how much space its    */
/* serialization will require.                                               */
/*                                                                           */
/* Arguments:                                                                */
/*    header to populate.                                                    */
/*    data container to populate.                                            */
/*    mxArray to analyze.                                                    */
/* Returns:                                                                  */
/*    size that shared memory segment will need to be.                       */
/* ------------------------------------------------------------------------- */
size_t deepscan(header_t *hdr, data_t *dat, const mxArray* mxInput) {
	/* counter */
	size_t i;

	/* for structures */	
	size_t field_num, count, plussiz;
	size_t strBytes;

	/* initialize data info; _possibly_ update these later */
	dat->pSize     = ( mwSize* )NULL;
	dat->pr        = (    void*)NULL;
	dat->pi        = (    void*)NULL;
	dat->ir        = (    void*)NULL;
	dat->jc        = ( mwIndex*)NULL;
	dat->child_dat = (  data_t*)NULL;
	dat->child_hdr = (header_t*)NULL;

	if ( mxInput==(const mxArray*)NULL || mxIsEmpty(mxInput) ) {
		/* initialize header info */
		hdr->isCell    = 0;
		hdr->isSparse  = 0;
		hdr->isComplex = 0;
		hdr->isStruct  = 0; 
		hdr->classid   = mxUNKNOWN_CLASS;
		hdr->nDims     = 0;
		hdr->elemsiz   = 0;
		hdr->nzmax     = 0;
		hdr->nFields   = 0; 
		hdr->shmsiz    = sizeof(header_t);
		return hdr->shmsiz;
	}

	/* initialize header info */
	hdr->isCell    = mxIsCell(mxInput);
	hdr->isSparse  = mxIsSparse(mxInput);
	hdr->isComplex = mxIsComplex(mxInput);
	hdr->isStruct  = mxIsStruct(mxInput);
	hdr->classid   = mxGetClassID(mxInput);
	hdr->nDims     = mxGetNumberOfDimensions(mxInput); 
	hdr->elemsiz   = mxGetElementSize(mxInput);
	hdr->nzmax     = mxGetNumberOfElements(mxInput);    /* update this later (if sparse) */
	hdr->nFields   = 0;									/* update this later */
	hdr->shmsiz    = sizeof(header_t);					/* update this later */

	/* copy the size */
	dat->pSize     = (mwSize*)mxGetDimensions(mxInput);	/* some abuse of const for now, fix on deep copy*/

	/* Add space for the dimensions */
	hdr->shmsiz += pad_to_align(hdr->nDims * sizeof(mwSize));

	if ( hdr->isStruct ) {
		/* struct case */

		mydebug("deepscan: found structure");

		/* How many fields to work with? */
		hdr->nFields = mxGetNumberOfFields(mxInput);

		/* Find the size required to store the field names */
		strBytes = FieldNamesSize(mxInput);     /* always a multiple of alignment size */
		hdr->shmsiz += strBytes;			    /* Add space for the field string */
		
		/* use mxCalloc so mem is freed on error via garbage collection */
		dat->child_hdr = (header_t*)mxCalloc(hdr->nzmax * hdr->nFields * (sizeof(header_t) + 
			sizeof(data_t)) + strBytes, 1); /* extra space for the field string */
		dat->child_dat = (data_t*)&dat->child_hdr[hdr->nFields * hdr->nzmax];
		dat->ir = (void*)&dat->child_dat[hdr->nFields * hdr->nzmax];

		/* make a record of the field names */
		CopyFieldNames(mxInput, (char*)(dat->ir));

		/* go through each recursively */
		count = 0;
		for ( i=0; i<hdr->nzmax; i++ ) {									/* each element */
			for (field_num = 0; field_num < hdr->nFields; field_num++) {	/* each field */
				/* call recursivley */
				plussiz      =	deepscan(&(dat->child_hdr[count]), 
			                    &(dat->child_dat[count]), mxGetFieldByNumber(mxInput,i,field_num) );
				hdr->shmsiz +=	plussiz;
				count++; /* progress */
				mydebug("deepscan: finished element %d field %d (%12u %12u)",i,field_num,plussiz,hdr->shmsiz);
			}
		}

	} else if ( hdr->isCell ) {
		/* cell case */

		mydebug("deepscan: found cell");

		/* use mxCalloc so mem is freed on error via garbage collection */
		dat->child_hdr = (header_t*)mxCalloc(hdr->nzmax,sizeof(header_t) + sizeof(data_t));
		dat->child_dat = (data_t*)&dat->child_hdr[hdr->nzmax];

		/* go through each recursively */
		for ( i=0;i<hdr->nzmax;i++ ) {
			plussiz      =	deepscan( &(dat->child_hdr[i]),
								&(dat->child_dat[i]), mxGetCell(mxInput,i));
			hdr->shmsiz +=	plussiz;
			mydebug("deepscan: finished %d (%12u %12u)",i,plussiz,hdr->shmsiz);
		}

		/* if its a cell it has to have at least one mom-empty element */
		if ( hdr->shmsiz==(1+hdr->nzmax)*sizeof(header_t) )
			mexErrMsgIdAndTxt("MATLAB:sharedmatrix:clone",
				"Required third argument (variable) must contain at "
				"least one non-empty item (at all levels).");

	} else if ( mxIsNumeric(mxInput) || mxIsLogical(mxInput) || mxIsChar(mxInput))  { /* a matrix containing data */
		/* matrix case */

		if ( hdr->isSparse ) {
			/* len(pr)=nzmax, len(ir)=nzmax, len(jc)=colc+1 */
			hdr->nzmax = mxGetNzmax(mxInput);/* make sure to do this because its sparse */
			dat->ir = mxGetIr(mxInput);
			dat->jc = mxGetJc(mxInput);
			/* ensure both pointers are aligned individually */
			hdr->shmsiz += pad_to_align(sizeof(mwIndex)*(hdr->nzmax));      /* ir */
			hdr->shmsiz += pad_to_align(sizeof(mwIndex)*(dat->pSize[1]+1)); /* jc */
		}

		dat->pr = mxGetData(mxInput);
		dat->pi = mxGetImagData(mxInput);
		/* ensure both pointers are aligned individually */
		hdr->shmsiz += pad_to_align(hdr->elemsiz*hdr->nzmax);		/* pr */
		if (hdr->isComplex)
			hdr->shmsiz += pad_to_align(hdr->elemsiz*hdr->nzmax);	/* pi */

	} else {
		/* error case */

		mexErrMsgIdAndTxt("MATLAB:sharedmatrix:clone",
			"Unknown class found: %s.", mxGetClassName(mxInput)); 
	}


	return hdr->shmsiz;
}


/* ------------------------------------------------------------------------- */
/* deepcopy                                                                  */
/*                                                                           */
/* Descend through header and data structure and copy relevent data to       */
/* shared memory.                                                            */
/*                                                                           */
/* Arguments:                                                                */
/*    header to serialize.                                                   */
/*    data container of pointers to data to serialize.                       */
/*    pointer to shared memory segment.                                      */
/* Returns:                                                                  */
/*    void                                                                   */
/* ------------------------------------------------------------------------- */
void deepcopy(header_t *hdr, data_t *dat, char *shm, bool allocate_only) {

	/* declarations */
	size_t i,offset;
	mwSize *pSize;			/* points to the size data */
	size_t nDims=hdr->nDims;

	/* for structures */	
	size_t nFields;
	int ret;

	/* copy the header (NOALIGN) */
	offset = sizeof(header_t);
	memcpy(shm,hdr,offset);
	shm+=offset;

	/* keep for later */
	pSize = (mwSize*)shm; 

	/* copy the dimensions (ALIGN) */
	offset = nDims * sizeof(mwSize);
	memcpy(shm,dat->pSize,offset);
	shm += pad_to_align(offset);

	if (hdr->isStruct) {
		/* struct case */

		mydebug("deepcopy: found structure");

		/* place the field names next in shared memory */
		/* note: NumFieldsFromString() automatically pads the offset */
		ret = NumFieldsFromString((char*)(dat->ir), &nFields, &offset);

		/* check the recovery */
		if ( ret || nFields!=hdr->nFields )
			mexErrMsgIdAndTxt("MATLAB:sharedmatrix:clone",
				"Structure fields have not been recovered properly.");

		/* copy the field string (ALIGN)  */
		memcpy(shm,dat->ir,offset);
		shm += offset; /* aligned already! in NumFieldsFromString */
		
		/* copy children recursively (NOALIGN) */
		for ( i=0; i<hdr->nzmax*hdr->nFields; i++ ) {
			deepcopy( &(dat->child_hdr[i]),&(dat->child_dat[i]), shm, allocate_only );
			shm += (dat->child_hdr[i]).shmsiz;
		}

	} else if ( hdr->isCell ) {
		/* cell case */ 

		mydebug("deepcopy: found cell");

		/* copy children recursively (NOALIGN) */
		for( i=0; i<hdr->nzmax; i++ ) {
			deepcopy(&(dat->child_hdr[i]),&(dat->child_dat[i]),shm, allocate_only );
			shm+=(dat->child_hdr[i]).shmsiz;
		}

	} else {
		/* matrix case */

		/* copy real data (ALIGN) */
		offset = (hdr->nzmax)*(hdr->elemsiz);
        if (!allocate_only)
            memcpy(shm,dat->pr,offset);
		shm += pad_to_align(offset);

		/* copy complex data (ALIGN) */
		if ( hdr->isComplex ) {
			offset = (hdr->nzmax)*(hdr->elemsiz);
           if (!allocate_only)
               memcpy(shm,dat->pi,offset);
			shm += pad_to_align(offset);
		}
		
		/* indices of sparse (ALIGN) */
		if ( hdr->isSparse ) {
			offset = hdr->nzmax*sizeof(mwIndex);
			memcpy(shm,dat->ir,offset);
			shm += pad_to_align(offset);

			offset = (pSize[1]+1)*sizeof(mwIndex);
			memcpy(shm,dat->jc,offset);
			shm += pad_to_align(offset);
		}
	}
}


/* ------------------------------------------------------------------------- */
/* deepfree                                                                  */
/*                                                                           */
/* Descend through header and data structure and free the memory.            */
/*                                                                           */
/* Arguments:                                                                */
/*    data container                                                         */
/* Returns:                                                                  */
/*    void                                                                   */
/* ------------------------------------------------------------------------- */
void deepfree(data_t *dat) {

	/* recurse to the bottom */
	if ( dat->child_dat!=(data_t*)NULL )
		deepfree(dat->child_dat);

	/* free on the way back up */
	mxFree(dat->child_hdr);

	dat->child_hdr = (header_t*)NULL;
	dat->child_dat = (data_t*)NULL;
	dat->ir        = (void*)NULL; /* could have been used for fieldNames */

}



/* ------------------------------------------------------------------------- */
/*                                                                           */
/* Andrew Smith Code for structs.                                            */
/*                                                                           */
/* ------------------------------------------------------------------------- */

/* Function to find the number of bytes required to store all of the */
/* field names of a structure									     */
int FieldNamesSize(const mxArray * mxStruct)
{
	const char*  pFieldName;
	int nFields;
	int i,j;			/* counters */
	int Bytes = 0;

	/* How many fields? */
	nFields = mxGetNumberOfFields(mxStruct);

	/* Go through them */
	for (i = 0; i < nFields; i++)
	{
		/* This field */
		pFieldName = mxGetFieldNameByNumber(mxStruct,i);
		j = 0;
		while (pFieldName[j])
		{	j++; }
		j++;		/* for the null termination */

		if (i == (nFields - 1))
			j++;	/* add this at the end for an end of string sequence since we use 0 elsewhere (term_char). */
		            

		/* Pad it out to 8 bytes */
		j = pad_to_align(j);

		/* keep the running tally */
		Bytes += j;

	}

	/* Add an extra alignment size to store the length of the string (in Bytes) */
	Bytes += align_size;

	return Bytes;

}

/* Function to copy all of the field names to a character array    */
/* Use FieldNamesSize() to allocate the required size of the array */
/* returns the number of bytes used in pList					   */
int CopyFieldNames(const mxArray * mxStruct, char* pList)
{
	const char*  pFieldName;    /* name of the field to add to the list */
	int nFields;				/* the number of fields */
	int i,j;					/* counters */
	int Bytes = 0;				/* number of bytes copied */

	/* How many fields? */
	nFields = mxGetNumberOfFields(mxStruct);

	/* Go through them */
	for (i = 0; i < nFields; i++)
	{
		/* This field */
		pFieldName = mxGetFieldNameByNumber(mxStruct,i);
		j = 0;
		while (pFieldName[j])
		{	pList[Bytes++]=pFieldName[j++];   }
		pList[Bytes++] = 0; /* add the null termination */

		/* if this is last entry indicate the end of the string */
		if (i == (nFields - 1))
		{	
			pList[Bytes++] = term_char;
			pList[Bytes++] = 0;  /* another null termination just to be safe */
		}

		/* Pad it out to 8 bytes */
		while (Bytes % align_size)
		{	pList[Bytes++] = 0;   }
	}

	/* Add an extra alignment size to store the length of the string (in Bytes) */
	*(unsigned int*)&pList[Bytes] = (unsigned int)(Bytes + align_size);
	Bytes += align_size;

	return Bytes;

}

/*Function to take point a each element of the char** array at a list of names contained in string */
/*ppCharArray must be allocated to length num_names */
/*names are seperated by null termination characters */
/*each name must start on an aligned address (see CopyFieldNames()) */
/*e.g. pCharArray[0] = &name_1, pCharArray[1] = &name_2 ...         */
/*returns 0 if successful */
int PointCharArrayAtString(char ** pCharArray, char* pString, int nFields)
{ 
	int i = 0;					    /* the index in pString we are up to */
	int field_num = 0;			    /* the field we are up to */

	/* and recover them */
	while (field_num < nFields)
	{
		/* scan past null termination chars */
		while (!pString[i])
			i++;

		/* Check the address is aligned (assume every forth bytes is aligned) */
		if (i % align_size)
			return -1;

		/* grab the address */
		pCharArray[field_num] = &pString[i];
		field_num++;

		/* continue until a null termination char */
		while (pString[i])
			i++;

	}
	return 0;
}

/* This function finds the number of fields contained within in a string */
/* the string is terminated by term_char, a character that can't be in a field name */
/* returns the number of field names found */
/* returns -1 if unaligned field names are found */
int NumFieldsFromString(const char* pString, size_t *pFields, size_t* pBytes)
{
	unsigned int stored_length;		/* the recorded length of the string */
	bool term_found = false;		/* has the termination character been found? */
	int i = 0;						/* counter */

	pFields[0] = 0;					/* the number of fields in the string */
	pBytes[0] = 0;


	/* and recover them */
	while (!term_found)
	{
		/* scan past null termination chars */
		while (!pString[i])
		{	i++;		 }

		/* is this the special termination character? */
		term_found = pString[i] == term_char;

		if (!term_found)
		{
			/* Check the address is aligned (assume every forth bytes is aligned) */
			if (i % align_size)
			{	return -1;    }

			/* a new field */
			pFields[0]++;

			/* continue until a null termination char */
			while (pString[i] && (pString[i] != term_char))
			{	i++;	}

			/* check there wasn't the terminal character immediately after the field */
			if (pString[i] == term_char)
			{	return -2;		}

		}
	}
	pBytes[0] = i;

	/* return the aligned size */
	pBytes[0] = pad_to_align(pBytes[0]);

	/* Check what's recovered matches the stored length */
	stored_length = *(unsigned int*)&pString[pBytes[0]] ;	  /* the recorded length of the string */
	pBytes[0] += align_size;								  /* add the extra space for the record */

	/* Check on it */
	if (stored_length != pBytes[0])
		return -3;

	return 0;

}

/* Function to find the bytes in the string starting from the end of the string (i.e. the last element is pString[-1])*/
/* returns < 0 on error */
int BytesFromStringEnd(const char* pString, size_t* pBytes)
{
	
	int offset = align_size;		  /* start from the end of the string */
	pBytes[0] = 0;					  /* default */

	/* Grab it directly from the string */
	pBytes[0] = (size_t)(*(unsigned int*)&pString[-offset]);

	/* scan for the termination character */
	while ((offset < 2*align_size) && (pString[-offset] != term_char))
	{ offset++;}

	if (offset == align_size)
	{	return -1;				}/* failure to find the control character */
	else
	{	return 0;				}/* happy */	 
	
}

