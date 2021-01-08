#ifndef TV_GPU_tex_HPP
#define TV_GPU_tex_HPP

#include "tmwtypes.h"


int checkLastError(char * msg);

void nonlocal_TV_init( float * Img, const float * Img0, const uint8_T *CloseInd,
       const uint8_T * neighbours, const float * p_Nbrs_weights,  const float dt, const float eps, const float lambda,const unsigned int Nclose, const unsigned int Nbrs, const int Rwin, 
       const unsigned int N, const unsigned int M, const unsigned int O, const unsigned int Niter);

void nonlocal_weight_TV_init( uint8_T *c, const float * p,
       const uint8_T * neighbours,  const unsigned int Nclose, const unsigned int Nclose_min,  const unsigned int Nbrs,  const int Rwin,  const int Rpatch, 
       const unsigned int N, const unsigned int M, const unsigned int O, const float threshold);

void local_TV_init( float * Img, const float dt, const float eps,  const unsigned int M, const unsigned int N, const unsigned int O, const unsigned int  Niter);

void local_TV_chambolle_init( float * Img, float ** Xi, const float dt, const float tau, 
                                const unsigned int M, const unsigned int N, const unsigned int O, const unsigned int  Niter);


#endif
