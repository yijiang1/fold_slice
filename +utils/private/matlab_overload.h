

#ifndef _MATLAB_OVERLOAD
#define _MATLAB_OVERLOAD

#include "mex.h"


// overload the GetData function for each of the possible data type + select the correct matlab get function 
inline void  GetData(const mxArray *in, const mxComplexDouble *& out) {  out =  mxGetComplexDoubles(in); return; }; 
inline void  GetData(const mxArray *in, const mxComplexSingle *& out) {  out =  mxGetComplexSingles(in); return;}; 
inline void  GetData(const mxArray *in, const mxComplexUint32 *& out) {  out =  mxGetComplexUint32s(in); return;}; 
inline void  GetData(const mxArray *in, const mxComplexUint16 *& out) {  out =  mxGetComplexUint16s(in); return;}; 
inline void  GetData(const mxArray *in, const mxComplexUint8 *& out)  {  out =  mxGetComplexUint8s(in); return;}; 
inline void  GetData(const mxArray *in, const mxDouble *& out)        {  out =  mxGetDoubles(in); return;}; 
inline void  GetData(const mxArray *in, const mxSingle *& out)        {  out =  mxGetSingles(in); return;}; 
inline void  GetData(const mxArray *in, const mxUint32 *& out)        {  out =  mxGetUint32s(in); return;}; 
inline void  GetData(const mxArray *in, const mxUint16 *& out)        {  out =  mxGetUint16s(in); return;}; 
inline void  GetData(const mxArray *in, const mxUint8 *& out)         {  out =  mxGetUint8s(in); return;}; 
inline void  GetData(const mxArray *in, mxComplexDouble *& out)       {  out =  mxGetComplexDoubles(in); return;}; 
inline void  GetData(const mxArray *in, mxComplexSingle *& out)       {  out =  mxGetComplexSingles(in); return;}; 
inline void  GetData(const mxArray *in, mxComplexUint32 *& out)       {  out =  mxGetComplexUint32s(in); return;}; 
inline void  GetData(const mxArray *in, mxComplexUint16 *& out)       {  out =  mxGetComplexUint16s(in); return;}; 
inline void  GetData(const mxArray *in, mxComplexUint8 *& out)        {  out =  mxGetComplexUint8s(in); return;}; 
inline void  GetData(const mxArray *in, mxDouble *& out)              {  out =  mxGetDoubles(in); return;}; 
inline void  GetData(const mxArray *in, mxSingle *& out)              {  out =  mxGetSingles(in); return;}; 
inline void  GetData(const mxArray *in, mxUint32 *& out)              {  out =  mxGetUint32s(in); return;}; 
inline void  GetData(const mxArray *in, mxUint16 *& out)              {  out =  mxGetUint16s(in); return;}; 
inline void  GetData(const mxArray *in, mxUint8 *& out)               {  out =  mxGetUint8s(in); return;}; 

inline void AddData_atomic( mxDouble &out, const mxDouble &in) {  
            #pragma omp atomic 
            out += in; }; 
inline void AddData_atomic( mxSingle &out, const mxSingle &in) {  
            #pragma omp atomic 
            out += in; }; 
inline void AddData_atomic( mxUint32 &out, const mxUint32 &in) {  
            #pragma omp atomic 
            out += in; }; 
inline void AddData_atomic( mxUint16 &out, const mxUint16 &in) {  
            #pragma omp atomic 
            out += in; }; 
inline void AddData_atomic( mxUint8  &out, const mxUint8  &in) {  
            #pragma omp atomic 
            out += in; }; 
inline void AddData_atomic( mxComplexDouble &out, const mxComplexDouble &in) { 
            #pragma omp atomic update
            out.real += in.real;
            #pragma omp atomic update
            out.imag += in.imag;}; 
inline void AddData_atomic( mxComplexSingle &out, const mxComplexSingle &in) { 
            #pragma omp atomic update
            out.real += in.real;
            #pragma omp atomic update
            out.imag += in.imag;}; 
inline void AddData_atomic( mxComplexUint32 &out, const mxComplexUint32 &in) { 
            #pragma omp atomic update
            out.real += in.real;
            #pragma omp atomic update
            out.imag += in.imag;}; 
inline void AddData_atomic( mxComplexUint16 &out, const mxComplexUint16 &in) {  
            #pragma omp atomic update
            out.real += in.real;
            #pragma omp atomic update
            out.imag += in.imag;}; 
inline void AddData_atomic( mxComplexUint8  &out, const mxComplexUint8 &in)  {  
            #pragma omp atomic update
            out.real += in.real;
            #pragma omp atomic update
            out.imag += in.imag;}; 

inline void AddData( mxDouble &out, const mxDouble &in) {  
            out += in; }; 
inline void AddData( mxSingle &out, const mxSingle &in) {  
            out += in; }; 
inline void AddData( mxUint32 &out, const mxUint32 &in) {  
            out += in; }; 
inline void AddData( mxUint16 &out, const mxUint16 &in) {  
            out += in; }; 
inline void AddData( mxUint8  &out, const mxUint8  &in) {  
            out += in; }; 
inline void AddData( mxComplexDouble &out, const mxComplexDouble &in) { 
            out.real += in.real; 
            out.imag += in.imag;}; 
inline void AddData( mxComplexSingle &out, const mxComplexSingle &in) { 
            out.real += in.real;
            out.imag += in.imag;}; 
inline void AddData( mxComplexUint32 &out, const mxComplexUint32 &in) { 
            out.real += in.real;
            out.imag += in.imag;}; 
inline void AddData( mxComplexUint16 &out, const mxComplexUint16 &in) {  
            out.real += in.real;
            out.imag += in.imag;}; 
inline void AddData( mxComplexUint8  &out, const mxComplexUint8 &in)  {  
            out.real += in.real;
            out.imag += in.imag;}; 
            
inline void SetData( mxDouble &out, const mxDouble &in) {  out = in; }; 
inline void SetData( mxSingle &out, const mxSingle &in) {  out = in; }; 
inline void SetData( mxUint32 &out, const mxUint32 &in) {  out = in; }; 
inline void SetData( mxUint16 &out, const mxUint16 &in) {  out = in; }; 
inline void SetData( mxUint8  &out, const mxUint8  &in) {  out = in; }; 
inline void SetData( mxComplexDouble &out, const mxComplexDouble &in) {  out.real = in.real; out.imag = in.imag; }; 
inline void SetData( mxComplexSingle &out, const mxComplexSingle &in) {  out.real = in.real; out.imag = in.imag; }; 
inline void SetData( mxComplexUint32 &out, const mxComplexUint32 &in) {  out.real = in.real; out.imag = in.imag; }; 
inline void SetData( mxComplexUint16 &out, const mxComplexUint16 &in) {  out.real = in.real; out.imag = in.imag; }; 
inline void SetData( mxComplexUint8  &out, const mxComplexUint8 &in)  {  out.real = in.real; out.imag = in.imag; }; 

#endif

