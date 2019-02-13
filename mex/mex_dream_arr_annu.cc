/***
*
* Copyright (C) 2003,2005,2006,2007,2008,2009,2014,2015 Fredrik Lingvall
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

// $Revision: 886 $ $Date: 2015-06-09 07:43:22 +0200 (Tue, 09 Jun 2015) $ $LastChangedBy: frli8848 $

#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <thread>
#include <mutex>
#include <signal.h>
#include "mex.h"
#include "dream_arr_annu.h"
#include "dream_error.h"

#define SINGLE 0
#define MULTIPLE 1

#ifdef USE_FFTW
#include "att.h"
#endif

//
// Globals
//
volatile int out_err = NONE;
std::mutex err_lock; 
int running;

//
// typedef:s
//

typedef struct
{
  size_t no;
  size_t start;
  size_t stop;
  double *RESTRICT ro;
  double dx;
  double dy;
  double dt; 
  size_t nt;
  int delay_method;
  double *RESTRICT delay;
  double v;
  double cp;
  double alfa;
  size_t isize;
  double *RESTRICT gr;
  int ifoc;
  int iweight;
  int iapo;
  double focal;
  double *RESTRICT apod;
  double param;
  double  *ud_focal;
  double *RESTRICT h; 
  int err_level;
} DATA;

typedef void (*sighandler_t)(int);

//
// Function prototypes.
//
void* smp_process(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

/***
 *
 * Thread function. 
 *
 ***/

void* smp_process(void *arg)
{
  int tmp_err = NONE, err = NONE;
  size_t n;
  DATA D = *(DATA *)arg;
  double xo, yo, zo;
  double *RESTRICT h = D.h;
  double dx=D.dx, dy=D.dy, dt=D.dt;
  size_t no=D.no, nt=D.nt;
  int    tmp_lev, err_level=D.err_level;
  double *RESTRICT delay=D.delay, *RESTRICT ro=D.ro, v=D.v, cp=D.cp, alfa=D.alfa;
  size_t start=D.start, stop=D.stop;
  int    ifoc=D.ifoc, iweight=D.iweight,iapo=D.iapo;
  double focal=D.focal, *apod=D.apod, param=D.param;
  size_t  isize = D.isize;
  double *RESTRICT gr=D.gr;
  double *RESTRICT ud_focal=D.ud_focal;

  // Let the thread finish and then catch the error. 
  if (err_level == STOP)
    tmp_lev = PARALLEL_STOP;
  else
    tmp_lev = err_level;

  if (ifoc != 6) { 
    
    if (D.delay_method == SINGLE) {
      for (n=start; n<stop; n++) {
	xo = ro[n]; 
	yo = ro[n+1*no]; 
	zo = ro[n+2*no]; 

	err = dream_arr_annu(xo,yo,zo,dx,dy,dt,nt,delay[0],v,cp,alfa,isize,gr,
			     ifoc,focal,apod,iweight,iapo,param,&h[n*nt],tmp_lev); 	

	if (err != NONE || out_err ==  PARALLEL_STOP) {
	  tmp_err = err;
	  if (err == PARALLEL_STOP || out_err ==  PARALLEL_STOP)
	    break; // Jump out when a STOP error occurs.
	}

	if (!running) {
	  std::cout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!" << std::endl;
	  return(NULL);
	}
	
      }
    } else { // MULTIPLE delays.
      for (n=start; n<stop; n++) {
	xo = ro[n]; 
	yo = ro[n+1*no]; 
	zo = ro[n+2*no]; 

	err = dream_arr_annu(xo,yo,zo,dx,dy,dt,nt,delay[n],v,cp,alfa,isize,
			     gr,ifoc,focal,apod,iweight,iapo,param,&h[n*nt],tmp_lev);
	
	if (err != NONE || out_err ==  PARALLEL_STOP) {
	  tmp_err = err;
	  if (err == PARALLEL_STOP || out_err ==  PARALLEL_STOP)
	    break; // Jump out when a STOP error occurs.
	}
	
	if (!running) {
	  std::cout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!" << std::endl;
	  return(NULL);
	}
	
      }
    }
  }
  else { // User defined focusing.
    
    if (D.delay_method == SINGLE) {
      for (n=start; n<stop; n++) {
	xo = ro[n]; 
	yo = ro[n+1*no]; 
	zo = ro[n+2*no]; 
	
	err = dream_arr_annu_ud(xo,yo,zo,dx,dy,dt,nt,delay[0],v,cp,alfa,isize,gr,
				ifoc,ud_focal,apod,iweight,iapo,param,&h[n*nt],tmp_lev); 	

	if (err != NONE || out_err ==  PARALLEL_STOP) {
	  tmp_err = err;
	  if (err == PARALLEL_STOP || out_err ==  PARALLEL_STOP)
	    break; // Jump out when a STOP error occurs.
	}

	if (!running) {
	  std::cout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!" << std::endl;
	  return(NULL);
	}

      }
    } else { // MULTIPLE delays.
      for (n=start; n<stop; n++) {
	xo = ro[n]; 
	yo = ro[n+1*no]; 
	zo = ro[n+2*no]; 
	
	err = dream_arr_annu_ud(xo,yo,zo,dx,dy,dt,nt,delay[n],v,cp,alfa,isize,
				gr,ifoc,ud_focal,apod,iweight,iapo,param,&h[n*nt],tmp_lev);
	
	if (err != NONE || out_err ==  PARALLEL_STOP) {
	  tmp_err = err;
	  if (err == PARALLEL_STOP || out_err ==  PARALLEL_STOP)
	    break; // Jump out when a STOP error occurs.
	}
	
	if (!running) {
	  std::cout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!" << std::endl;
	  return(NULL);
	}

      }
    }
  }

  // Lock out_err for update, update it, and unlock. 
  err_lock.lock();
  
  if ((tmp_err != NONE) && (out_err == NONE))
    out_err = tmp_err;
  
  err_lock.unlock();
  
  return(NULL);
}

/***
 *
 * Signal handlers.
 *
 ***/

void sighandler(int signum) {
  //printf("Caught signal SIGTERM.\n");
  running = false;
}

void sig_abrt_handler(int signum) {
  //printf("Caught signal SIGABRT.\n");
}

void sig_keyint_handler(int signum) {
  //printf("Caught signal SIGINT.\n");
}


/***
 *
 * Matlab (MEX) gateway function for dream_arr_annu
 *
 ***/

extern void _main();

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *RESTRICT ro, *RESTRICT s_par,*RESTRICT m_par;
  char   apod_met[50],foc_met[50];
  int    buflen;
  double dx,dy,dt;
  size_t nt, no;
  double param=0,*RESTRICT delay,v,cp,alfa;
  size_t isize=0;
  double *RESTRICT gr;
  int    ifoc=0;
  double focal=0,*RESTRICT ud_focal=NULL;
  double *RESTRICT apod=NULL;
  int    iweight=0, iapo=0;
  double *RESTRICT h, *err_p;
  int    err_level=STOP, is_set = false;
  char   err_str[50];
  DATA   *D;
  size_t start, stop;
  std::thread *threads;
  unsigned int thread_n, nthreads;
  sighandler_t old_handler, old_handler_abrt, old_handler_keyint;
  
  // Check for proper number of arguments.
  if (!((nrhs == 10) || (nrhs == 11))) {
    dream_err_msg("dream_arr_annu requires 10 or 11 input arguments!");
  }
  else
    if (nlhs > 2) {
      dream_err_msg("Too many output arguments for dream_arr_annu!");
    }

  //
  // Observation point.
  //
  
  // Check that arg (number of observation points) x 3 matrix
  if (mxGetN(prhs[0]) != 3)
    dream_err_msg("Argument 1 to dream_arr_annu must be a (number of observation points) x 3 matrix!");
  
  no = mxGetM(prhs[0]); // Number of observation points.
  ro = mxGetPr(prhs[0]);

  //
  // Grid function (position vectors of the elements).
  //
  
  if ((mxGetM(prhs[1]) > 1) & (mxGetN(prhs[1]) > 1))
    dream_err_msg("Argument 2 must a vector (number of array elements)");
  
  isize = (int) mxGetM(prhs[1])*mxGetN(prhs[1]); // Number of elementents in the array.
  gr = mxGetPr(prhs[1]);	// Vector of annular radi,
  
  //
  // Temporal and spatial sampling parameters.
  //
  
 
  // Check that arg 2 is a 4 element vector
  if (!((mxGetM(prhs[2])==4 && mxGetN(prhs[2])==1) || (mxGetM(prhs[2])==1 && mxGetN(prhs[2])==4)))
    dream_err_msg("Argument 3 must be a vector of length 4!");
  
  s_par = mxGetPr(prhs[2]);
  dx    = s_par[0];		// Spatial x discretization size.
  dy    = s_par[1];		// Spatial dy iscretization size.
  dt    = s_par[2];		// Temporal discretization size (= 1/sampling freq).
  nt    = (size_t) s_par[3];	// Length of SIR.
  
  //
  // Start point of impulse response vector ([us]).
  //
  
  // Check that arg 4 is a scalar.
  if ( (mxGetM(prhs[3]) * mxGetN(prhs[3]) !=1) && ((mxGetM(prhs[3]) * mxGetN(prhs[3])) != no))
    dream_err_msg("Argument 4 must be a scalar or a vector with a length equal to the number of observation points!");
  
  delay = mxGetPr(prhs[3]);
  
  //
  // Material parameters
  //
  
  // Check that arg 5 is a 3 element vector.
  if (!((mxGetM(prhs[4])==3 && mxGetN(prhs[4])==1) || (mxGetM(prhs[4])==1 && mxGetN(prhs[4])==3)))
    dream_err_msg("Argument 5 must be a vector of length 3!");
  
  m_par = mxGetPr(prhs[4]);
  v     = m_par[0]; // Normal velocity of transducer surface.
  cp    = m_par[1]; // Sound speed.
  alfa  = m_par[2]; // Attenuation coefficient [dB/(cm MHz)],

  //  
  // Focusing parameters.
  //

  //  ifoc = 1 - no foc, 2 foc xy.

  if (nrhs >= 6) {
 
    if (!mxIsChar(prhs[5]))
      dream_err_msg("Argument 6 must be a string");
    
    buflen = (mxGetM(prhs[5]) * mxGetN(prhs[5]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[5],foc_met,buflen);
    
    if (!strcmp(foc_met,"off"))
      ifoc = 1; 
    
    if (!strcmp(foc_met,"on"))
      ifoc = 2; 
    
    if (!strcmp(foc_met,"ud")) {
      ifoc = 6; 
      is_set = true;
      
      if (mxGetM(prhs[6]) * mxGetN(prhs[6]) != isize ) {
	printf("The time delay vector (argument 7) for user defined ('ud') focusing\n") ;
	dream_err_msg("delays must have the same length as the number of array elements.!");
	
      }
      ud_focal = mxGetPr(prhs[6]);
    }
    else {
      
      // Check that arg 7 is a scalar.
      if (mxGetM(prhs[6]) * mxGetN(prhs[6]) !=1 )
	dream_err_msg("Argument 7 must be a scalar!");
      
      // Focal point (in mm).
      focal = mxGetScalar(prhs[6]);
    }

  } else
    ifoc = 1;
  
  //
  // Apodization.
  //

  // iweight = 1 - no weighting, 2  weighting.
  // iapo = 0 - user defined, 1 traingle, 2 Gauss, 3 raised cosine, 4 simply supported, 5 clamped. 

  if (nrhs >= 8) {
    
    if (!mxIsChar(prhs[7]))
      dream_err_msg("Argument 8 must be a string");
    
    
    buflen = (mxGetM(prhs[7]) * mxGetN(prhs[7]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[7],apod_met,buflen);
    
    iweight = 1;			// default off. 
    is_set = false;
    
    if (!strcmp(apod_met,"off")) {
      iweight = 1;
      is_set = true; 
    }
    
    if (!strcmp(apod_met,"ud")) {
      iweight = 2; 
      iapo = 0; 
      is_set = true; 

      // Vector of apodization weights.
      if (mxGetM(prhs[8]) * mxGetN(prhs[8]) != isize)
	dream_err_msg("The length of argument 9 (apodization vector) must be the same as the number of array elements!");
      
      apod = mxGetPr(prhs[8]); 
    }
    
    if (!strcmp(apod_met,"triangle")) {
      iweight = 2; 
      iapo = 1; 
      is_set = true; 
    }
    
    if (!strcmp(apod_met,"gauss")) {
      iweight = 2; 
      iapo = 2; 
      is_set = true; 
    }
    
    if (!strcmp(apod_met,"raised")) {
      iweight = 2; 
      iapo = 3; 
      is_set = true; 
    }
    
    if (!strcmp(apod_met,"simply")) {
      iweight = 2; 
      iapo = 4; 
      is_set = true; 
    }
    
    if (!strcmp(apod_met,"clamped")) {
      iweight = 2; 
      iapo = 5; 
      is_set = true; 
    }

    if (is_set == false)
      dream_err_msg("Unknown apodization method!");
    
    // Parameter for raised cos and Gaussian apodization functions.
    if (mxGetM(prhs[9]) * mxGetN(prhs[9]) !=1 )
      dream_err_msg("Argument 10 must be a scalar");
    
    param = mxGetScalar(prhs[9]);  
  }
  else 
    iweight = 1;   

  //
  // Number of threads.
  //

  // Get number of CPU cores (including hypethreading, C++11)
  nthreads = std::thread::hardware_concurrency();
  
  // nthreads can't be larger then the number of observation points.
  if (nthreads > no) { 
    nthreads = no;
  }

  //
  // Error Reporting.
  //

  if (nrhs == 11) {

    if (!mxIsChar(prhs[10]))
      dream_err_msg("Argument 11 must be a string");
    
    buflen = (mxGetM(prhs[10]) * mxGetN(prhs[10]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[10],err_str,buflen);

    is_set = false;
    
    if (!strcmp(err_str,"ignore")) {
      err_level = IGNORE; 
      is_set = true;
    }
    
    if (!strcmp(err_str,"warn")) {
      err_level = WARN; 
      is_set = true;
    }

    if (!strcmp(err_str,"stop")) {
      err_level = STOP; 
      is_set = true;
    }

    if (is_set == false)
      dream_err_msg("Unknown error level!");
  }
  else
    err_level = STOP; // Default.
  
  // Create an output matrix for the impulse response
  plhs[0] = mxCreateDoubleMatrix(nt,no,mxREAL);
  h = mxGetPr(plhs[0]);

  //
  // Register signal handlers.
  //

  if ((old_handler = signal(SIGTERM, &sighandler)) == SIG_ERR) {
    printf("Couldn't register SIGTERM signal handler.\n");
  }

  if (( old_handler_abrt=signal(SIGABRT, &sighandler)) == SIG_ERR) {
    printf("Couldn't register SIGABRT signal handler.\n");
  }
  
  if (( old_handler_keyint=signal(SIGINT, &sighandler)) == SIG_ERR) {
    printf("Couldn't register SIGINT signal handler.\n");
  }
  
  //
  // Call the DREAM subroutine.
  //

  out_err = NONE;
  running = true;

#ifdef USE_FFTW
  if (alfa != (double) 0.0)
    att_init(nt,nthreads);
#endif

  // Allocate local data.
  D = (DATA*) malloc(nthreads*sizeof(DATA));

  // Allocate mem for the threads.
  threads = new std::thread[nthreads]; // Init thread data.
  
  for (thread_n = 0; thread_n < nthreads; thread_n++) {

    start = thread_n * no/nthreads;
    stop =  (thread_n+1) * no/nthreads;

    // Init local data.
    D[thread_n].start = start; // Local start index;
    D[thread_n].stop = stop; // Local stop index;
    D[thread_n].no = no;
    D[thread_n].ro = ro;
    D[thread_n].dx = dx;
    D[thread_n].dy = dy;
    D[thread_n].dt = dt; 
    D[thread_n].nt = nt;

    if (mxGetM(prhs[3]) * mxGetN(prhs[3]) == 1)
      D[thread_n].delay_method = SINGLE; // delay is a scalar.
    else
      D[thread_n].delay_method = MULTIPLE; // delay is a vector.

    D[thread_n].delay = delay;
    D[thread_n].v = v;
    D[thread_n].cp = cp;
    D[thread_n].alfa = alfa;
    D[thread_n].isize = isize;
    D[thread_n].gr = gr;
    D[thread_n].ifoc = ifoc;
    D[thread_n].iweight = iweight;
    D[thread_n].iapo = iapo;
    D[thread_n].focal = focal;
    D[thread_n].apod = apod;
    D[thread_n].param = param;
    D[thread_n].ud_focal = ud_focal;
    D[thread_n].h = h;
    D[thread_n].err_level = err_level;
    
    // Starts the threads.
    threads[thread_n] = std::thread(smp_process, &D[thread_n]); // Start the threads.
  }
  
  // Wait for all threads to finish.
  for (thread_n = 0; thread_n < nthreads; thread_n++)
    threads[thread_n].join();
  
  // Free memory.
  free((void*) D);
 
  //
  // Restore old signal handlers.
  //

  if (signal(SIGTERM, old_handler) == SIG_ERR) {
    printf("Couldn't register old SIGTERM signal handler.\n");
  }
   
  if (signal(SIGABRT,  old_handler_abrt) == SIG_ERR) {
    printf("Couldn't register old SIGABRT signal handler.\n");
  }
  
  if (signal(SIGINT, old_handler_keyint) == SIG_ERR) {
    printf("Couldn't register old SIGINT signal handler.\n");
  }

#ifdef USE_FFTW
  if (alfa != (double) 0.0)
    att_close();
#endif 

  if (!running) {
    dream_err_msg("CTRL-C pressed!\n"); // Bail out.
  }

  //
  // Check for Error.
  //

  if ( (err_level == STOP) && (out_err != NONE))
    dream_err_msg(""); // Bail out if error.
 
  //  
  // Return error.
  //
  if (nlhs == 2) {
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    err_p =  mxGetPr(plhs[1]);
    err_p[0] = (double) out_err;
  }
  
  return;
}
      
