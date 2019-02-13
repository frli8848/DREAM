
/***
*
* Copyright (C) 2010,2012,2014,2015 Fredrik Lingvall
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

// $Revision: 857 $ $Date: 2015-05-05 21:22:14 +0200 (Tue, 05 May 2015) $ $LastChangedBy: frli8848 $

#include <thread> // C++11 threads.

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "affinity.h"

#ifdef __gnu_linux__
//#define _GNU_SOURCE // Note this is a nonstandard GNU extensions.
#include <pthread.h>
#include <sched.h>
#endif

/***
 *
 * The affinity can be set with the taskset
 * command as well.
 *
 * fredrik@smodcomp:~> taskset -p 1 11623
 * pid 11623's current affinity mask: 2
 * pid 11623's new affinity mask: 1
 *
 ***/

int set_dream_thread_affinity(dream_idx_type thread_n, dream_idx_type nthreads, std::thread *threads)
{
  int err = 0;
  
#ifdef __gnu_linux__   

  // Linux uses pthreads.
  pthread_t my_threads = threads->native_handle(); // Get the (implementation defined) underlying thread handle
  
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset); // Clear the cpu set.

  dream_idx_type cpu;

  // CPU affinity for the n:th thread.
  err = pthread_getaffinity_np(my_threads, sizeof(cpu_set_t),&cpuset);
  if (err != 0) {
    perror("Error in affinity.c when reading cpu affinity\n");
    return err;
  }

#ifdef DEBUG
  printf("\nInitial CPU affinity mask for thread is : | ");    
  for (cpu = 0; cpu < nthreads; cpu++)
    printf("%d | ", CPU_ISSET(cpu, &cpuset));
  printf("\n");
#endif

  //    
  // Try to set/clear the affinity mask to allow/disallow to run on a CPU/Core.
  //

  // Allow thread number 'thread_n' to run on core (cpu) number 'thread_n' only.
  for (cpu=0; cpu < nthreads; cpu++) {
    if (cpu == thread_n)
      CPU_SET(cpu, &cpuset);
    else
      CPU_CLR(cpu, &cpuset);
  }
  
  // Clear the remaining cpu set.
  for (cpu=nthreads; cpu < CPU_SETSIZE; cpu++) 
    CPU_CLR(cpu, &cpuset);
  
  err = pthread_setaffinity_np(my_threads, sizeof(cpu_set_t),&cpuset);
  if (err != 0)
    perror("Error in affinity.c when setting cpu affinity\n");
  
#ifdef DEBUG
  printf("\n        CPU affinity mask for thread is : | ");    
  for (cpu = 0; cpu < nthreads; cpu++)
    printf("%d | ", CPU_ISSET(cpu, &cpuset));
  printf("\n");
#endif

#endif
  
  return err;
}
