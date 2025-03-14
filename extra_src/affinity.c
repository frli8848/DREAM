/***
*
* Copyright (C) 2010,2012,2014 Fredrik Lingvall
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


// Note this is a nonstandard GNU extensions.
#define _GNU_SOURCE
#include <pthread.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

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

int set_dream_thread_affinity(int max_cpu_num, pthread_t *threads, int *affinity_mask)
{
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset); // Clear the cpu set.
  int err = 0;
  size_t cpu;

  // CPU affinity for the n:th thread.
  err = pthread_getaffinity_np(*threads, sizeof(cpu_set_t),&cpuset);
  if (err != 0) {
    perror("Error in affinity.c when reading cpu affinity\n");
    return err;
  }

#ifdef DEBUG
  printf("\nInitial CPU affinity mask for thread is : | ");
  for (cpu = 0; cpu < max_cpu_num; cpu++)
    printf("%d | ", CPU_ISSET(cpu, &cpuset));
  printf("\n");
#endif

  //
  // Try to set/clear the affinity mask to allow/disallow to run on a CPU/Core.
  //

  // Allow only those CPUs set in the affinity_mask.
  for (cpu=0; cpu < max_cpu_num; cpu++) {
    if (affinity_mask[cpu] > 0)
      CPU_SET(cpu, &cpuset);
    else
      CPU_CLR(cpu, &cpuset);
  }

  // Clear the remaining cpu set.
  for (cpu=max_cpu_num; cpu < CPU_SETSIZE; cpu++)
    CPU_CLR(cpu, &cpuset);

  err = pthread_setaffinity_np(*threads, sizeof(cpu_set_t),&cpuset);
  if (err != 0)
    perror("Error in affinity.c when setting cpu affinity\n");

#ifdef DEBUG
  printf("\n        CPU affinity mask for thread is : | ");
  for (cpu = 0; cpu < max_cpu_num; cpu++)
    printf("%d | ", CPU_ISSET(cpu, &cpuset));
  printf("\n");
#endif

  return err;
}
