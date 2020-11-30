/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * 
 * @section LICENSE
 * Copyright (c) 2019, SeisSol Group
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * 
 **/

#include "Pin.h"

#include <iostream>
#include <sys/sysinfo.h>
#include <sched.h>
#include <sstream>
#include <omp.h>

seissol::parallel::Pinning::Pinning() {
  CPU_ZERO(&pthreadMask);

  // Affinity mask of the entire process
  sched_getaffinity(0, sizeof(cpu_set_t), &processMask);
}
void seissol::parallel::Pinning::init() {
  // Affinity mask for the OpenMP workers
  openmpMask = getWorkerUnionMask();

  // Affinity mask of the pthreads -> all free cores
  CPU_XOR(&pthreadMask, &processMask, &openmpMask);
}

cpu_set_t seissol::parallel::Pinning::getWorkerUnionMask() const {
  cpu_set_t workerUnion;
  CPU_ZERO(&workerUnion);
#ifdef _OPENMP
  #pragma omp parallel
  {
    cpu_set_t worker;
    CPU_ZERO(&worker);
    sched_getaffinity(0, sizeof(cpu_set_t), &worker);
    #pragma omp critical
    {
      CPU_OR(&workerUnion, &workerUnion, &worker);
    }
  }
#else
  sched_getaffinity(0, sizeof(cpu_set_t), &workerUnion);
#endif

  return workerUnion;
}

cpu_set_t seissol::parallel::Pinning::getFreeCPUsMask() const {
  return pthreadMask;
}

bool seissol::parallel::Pinning::freeCPUsMaskEmpty(cpu_set_t const& set) {
  return CPU_COUNT(&set) == 0;
}

void seissol::parallel::Pinning::pinToFreeCPUs() const {
  sched_setaffinity(0, sizeof(cpu_set_t), &pthreadMask);
}

std::string seissol::parallel::Pinning::maskToString(cpu_set_t const& set) {
  std::stringstream st;
  for (int cpu = 0; cpu < get_nprocs(); ++cpu) {
    if (cpu % 10 == 0 && cpu != 0 && cpu != get_nprocs()-1) {
      st << '|';
    }
    if (CPU_ISSET(cpu, &set)) {
      st << cpu % 10;
    } else {
      st << '-';
    }
  }
  return st.str();
}
