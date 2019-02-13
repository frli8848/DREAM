/***
*
* Copyright (C) 2005,2006,2007,2008 Fredrik Lingvall
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

//#include "dreamrect.h"
//#include "dream_error.h"

// Error levels.
#define NONE   0
#define STOP   -1
#define WARN   -2
#define IGNORE -3
#define PARALLEL_STOP -4

#define TRUE 1
#define FALSE 0

/***
 *
 * test_dreamrect.c -  for evaluating compiler flags for DREAM.
 *
 ***/


int dream_out_of_bounds_err(char *msg, int err, int err_level);
void dream_err_msg(char *msg);


#define mexErrMsgTxt printf  


int dream_out_of_bounds_err(char *msg, int err, int err_level)
{

  switch(err_level) {

  case STOP:
    printf(msg);
    printf(" (offset = %d samples)\n",err);
    //mexErrMsgTxt(""); // Bail out!
    exit(0);
    break;

  case WARN:
    printf(msg);
    printf(" (offset = %d samples)\n",err);
    break;

  case IGNORE:
    break; // Do nothing.

  case PARALLEL_STOP:
    printf(msg);
    printf(" (offset = %d samples)\n",err);
    break; 

  default:
    break;
  }
    
  return err_level;
}

void dream_err_msg(char *msg)
{
  mexErrMsgTxt(msg);

  return;
}


//
// Dummy functions.
//

/***
 *
 * Amplitude spectrum of a real vector.
 *
 ***/

void dream_fft(double *x, double *y, int n)
{
  return;
}

/***
 *
 * Inverse DFT of a real vector (returns the real part).
 *
 ***/

void dream_ifft(double *x, double *y, int n)
{

  return;
}

/***
 *
 * Inverse DFT of a complex vector (returns the real part).
 *
 ***/

void cr_ifft(double *xir, double *xii, double *y, int n)
{

  return;
}
/***
 *
 * Inverse DFT of a complex vector.
 *
 ***/

void cc_ifft(double *xir, double *xii, double *yor, double *yoi, int n)
{

  return;
}


/***
 *
 * dreamrect -
 *
 ***/

int dreamrect(double xo,
	      double yo,
	      double zo,
	      double a,
	      double b,
	      double dx,
	      double dy,
	      double dt,
	      int nt,
	      double delay,
	      double v,
	      double cp,
	      double alfa,
	      double *h,
	      int err_level);


//
//  Function prototypes.
//

void modri(double xo,
	   double yo,
	   double zo,
	   double x,
	   double y,
	   double *ri);

/***
 *
 * subroutine dreamrect
 *
 ***/
 
int dreamrect(double xo,
	       double yo,
	       double zo,
	       double a,
	       double b,
	       double dx,
	       double dy,
	       double dt,
	       int nt,
	       double delay,
	       double v,
	       double cp,
	       double alfa,
	       double *h,
	       int err_level)
{
  int i, j, it;
  double t;
  double ai, ds, pi;
  double ri, tt, qan, xsi, ysj;
  double xsmin, xsmax, ysmin, ysmax;
  int err = NONE;

  xsmin = -a/2;  
  xsmax = a/2;
  ysmin = -b/2;
  ysmax = b/2;

  pi = atan( (double) 1.0) * 4.0;
  ds = dx * dy;

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0 ;
  }
  
  j = 0;
  j++;
  ysj = ysmin + (j-1) * dy + dy / 2;
  while (ysj <= ysmax) {
    
    i = 0;  
    i++;
    xsi = xsmin + (i-1) * dx + dx / 2;
    while (xsi <= xsmax) {
      
      modri(xo, yo, zo, xsi, ysj, &ri);
      ai = v * ds / (2*pi * ri);
      ai /= dt;
      /* en SI */
      ai *= 1000;
      /* del prop en micros */
      t = ri * 1000/cp;
      tt = t - delay;
      qan = tt / dt;
      it = (int) rint(qan);

      // Check if index is out of bounds.
      if ( (it < nt) && (it >= 0) ) {
	
	// Check if absorbtion is present.
	if (alfa == (double) 0.0) {
	  h[it] += ai;
	} else {;
	//att(alfa,ri,it,dt,cp,h,nt,ai); Not used here.
	}
      }
      else {
	if  (it >= 0)
	  err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
	else
	  err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);

	if (err_level == PARALLEL_STOP)
	  return err; // Bail out.
      }
      
      i++;
      xsi = xsmin + (i-1)*dx + dx/2;
    } 
    
    j++;
    ysj = ysmin + (j-1)*dy + dy / 2;
  }
  
  return err;
} /* dreamrect */


/***
 * 
 * subrutine modri(xi,xs,hs,ri,rx,rz) pour trouver le longeur du vecteur 
 *
 ***/

void modri(double xo,
	   double yo,
	   double zo,
	   double x,
	   double y,
	   double *ri)
{
  double rx, ry, rz;

  ry = yo - y;
  rx = xo - x;
  rz = zo;
  *ri = sqrt(rx*rx + rz*rz + ry*ry);
  
  return;
} /* modri */


int main(int argc, char ** argv)
{
  int    nt=0,no,n;
  double xo,yo,zo, a, b, dx, dy, dt;
  double delay,v,cp,alfa;     
  double *h,sum;
  int    err_level=STOP, err=NONE;

  int i, number_of_SIRs;

  number_of_SIRs = 1000;	// Takes in the order of one second (a = 10.0,
				// b = 15.0).


  // do we have verbose output?
  bool ga_testing = false;
  
  if (argc > 1) {
    for (i = 1; i < argc; ++i) {
      if (!strcmp(argv[1],"-ga")) {
	ga_testing = true;
	break;
      }
    }
  }
  
  
  cp = 1000.0;

  a = 10.0;
  b = 15.0;
  //b = 15.0/2;

  nt = 2500;
  //nt = 25000; // Set this to a high value to generate cache misses.
  dx = 0.1;
  dy = 0.1;
  dt = 1.0/40.0;

  alfa = 0;
  v = 1.0;
  delay = 10.0/ cp * 1000.0;


  // Create output space for the impulse response(s).
  h = (double*) malloc(number_of_SIRs * nt * sizeof(double));
  if (!h) {
    printf("Memory allocation failed!\n");
    exit(0);
  }


  // Start timing.    
  struct timespec start, stop;
  clock_gettime(CLOCK_REALTIME,&start);
  
  // Call the DREAM subroutine.
  no = number_of_SIRs;
  //no = 1000;

  for (n=0; n<no; n++) {
    xo = ((double) n) / 100.0;
    yo = 0.0; 
    zo = 10.0;
    //err = dreamrect(xo,yo,zo,a,b,dx,dy,dt,nt,delay,v,cp,alfa,h,err_level);    
    err = dreamrect(xo,yo,zo,a,b,dx,dy,dt,nt,delay,v,cp,alfa,&h[n*nt],err_level);    
  }

  // Calculate run time.
  clock_gettime(CLOCK_REALTIME,&stop);        
  double run_time = (stop.tv_sec - start.tv_sec) + (double)(stop.tv_nsec - start.tv_nsec) / 1000000000.0;


  
  // Report runtime (to Acovea).
  if (ga_testing)
    fprintf(stdout,"%f",run_time);
  else {
    sum = 0.0;
    
    for (n=0; n<(nt*number_of_SIRs);n++)
      sum += h[n];

    fprintf(stdout,"\ndreamrect run time: %f [s] (sum = %2.10f)\n\n",run_time);
  }


  free( (void*) h);

  fflush(stdout);
    
  return(0);
}
      

