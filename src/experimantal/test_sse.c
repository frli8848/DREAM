/***
*
* Copyright (C) 2006 Fredrik Lingvall
*
* This file is part of the DREAM Toolbox.
*
* The DREAM Toolbox is free software; you can redistribute it and/or modify 
* it under the terms of the GNU General Public License as published by the
* Free Software Foundation; either version 2, or (at your option) any
* later version.
*
* The DREAM Toolbox is distributed in the hope that it will be useful, but 
* WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
* for more details.
*
* You should have received a copy of the GNU General Public License
* along with the DREAM Toolbox; see the file COPYING.  If not, write to the 
* Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
* 02110-1301, USA.
*
***/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>

#include <stdbool.h>
#include <time.h>
#include <stdint.h>

// SSE
#ifdef  __SSE2__
//#include <pmmintrin.h>
//#include <xmmintrin.h>
#include <emmintrin.h>

//typedef double v2df __attribute__ ((mode(V2DF)));

//typedef double double_sse  __attribute__((aligned(16)));

//#define double_sse v2df

#define AlignData(data) (void *)(( (int) data + 15) &~ 0x0F)

#endif


int main() 
{
  double *a, *b, *c1, *c2, *c22 , *a2, *b2;
  double err = 0.0;
  int n,N;
  struct timespec start, stop;

  N = 2*10000000;
  //N = 100;

  a = (double*) malloc(N * sizeof(double));
  b = (double*) malloc(N * sizeof(double));
  c1 = (double*) malloc(N * sizeof(double));

  for (n=0; n<N; n++) {
    a[n] = (double) n * 1.2;
    b[n] = (double) -n * 2;
  }


  printf("%p\n",a);


  printf("Usual way...\n");

  // Start timing.    
  clock_gettime(CLOCK_REALTIME,&start);

  for (n=0; n<N; n++) {
    c1[n] = a[n] * b[n];
  }

  // Calculate run time.
  clock_gettime(CLOCK_REALTIME,&stop);        
  double run_time = (stop.tv_sec - start.tv_sec) + (double)(stop.tv_nsec - start.tv_nsec) / 1000000000.0;
  fprintf(stdout,"run time: %f\n\n",run_time);
  

#ifdef __SSE__

  printf("...with SSE2 intrinsics...\n");
  
  __m128d a_sse, b_sse, c_sse;

  c2 = (double*) malloc(N * sizeof(double));

  a2= AlignData(a);
  a_sse = _mm_load_pd(a2);
  b2= AlignData(b);
  b_sse = _mm_load_pd(b2);
  c22= AlignData(c2);
  c_sse = _mm_load_pd(c22);

  clock_gettime(CLOCK_REALTIME,&start);



  // For non 16-bit aligned data.
  //for (n=0; n<N-1; n+=2) {
  //  a_sse = _mm_loadu_pd(&a[n]);
  //  b_sse = _mm_loadu_pd(&b[n]);
  //  c_sse = _mm_mul_pd(a_sse,b_sse);
  //  _mm_storeu_pd(&c2[n],c_sse);
  //}

  // For 16-bit aligned data.
  for (n=0; n<N-1; n+=2) {
    a_sse = _mm_load_pd(&a2[n]);
    b_sse = _mm_load_pd(&b2[n]);
    c_sse = _mm_mul_pd(a_sse,b_sse);
    _mm_store_pd(&c22[n],c_sse);
  }

  // Calculate run time.
  clock_gettime(CLOCK_REALTIME,&stop);        
  run_time = (stop.tv_sec - start.tv_sec) + (double)(stop.tv_nsec - start.tv_nsec) / 1000000000.0;
  fprintf(stdout,"run time: %f\n\n",run_time);

  for (n=0; n<N; n++) {
    err += sqrt(c1[n]*c1[n] - c2[n]*c2[n]);
  }

  fprintf(stdout,"err =  %f\n\n",err);

  free( (void*) c2);

#endif


  free( (void*) a);
  free( (void*) b);
  free( (void*) c1);
    
  return 0;
}
