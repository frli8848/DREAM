/***
*
* Copyright (C) 2010,2012,2014,2015,2019,2021 Fredrik Lingvall
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

#include <thread> // C++11 threads.

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "affinity.h"

// NB. We only support pthreads for setting the affinity

#ifdef HAVE_PTHREADS
#include <pthread.h>
#include <sched.h>
#endif

// NB. The affinity can be set with the taskset command as well.

int set_dream_thread_affinity(dream_idx_type thread_n, dream_idx_type nthreads, std::thread *threads)
{
  int err = 0;

#ifdef HAVE_PTHREADS

  cpu_set_t cpuset;

#ifdef DEBUG
  printf("\nInitial CPU affinity mask for thread is : | ");
  for (cpu = 0; cpu < nthreads; cpu++)
    printf("%d | ", CPU_ISSET(cpu, &cpuset));
  printf("\n");
#endif

  CPU_ZERO(&cpuset);
  CPU_SET(thread_n, &cpuset);
  err = pthread_setaffinity_np(threads[thread_n].native_handle(), sizeof(cpu_set_t), &cpuset);

#endif

  if (err != 0) {
    perror("Error in set_dream_thread_affinity when setting the CPU affinity\n");
    return err;
  }

  return err;
}
