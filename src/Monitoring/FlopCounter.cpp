/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013, SeisSol Group
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
 * Counts the floating point operations in SeisSol.
 **/

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "FlopCounter.hpp"

#include <utils/logger.h>

// Define the FLOP counter.
long long libxsmm_num_total_flops = 0;

long long g_SeisSolNonZeroFlopsLocal = 0;
long long g_SeisSolHardwareFlopsLocal = 0;
long long g_SeisSolNonZeroFlopsNeighbor = 0;
long long g_SeisSolHardwareFlopsNeighbor = 0;
long long g_SeisSolNonZeroFlopsOther = 0;
long long g_SeisSolHardwareFlopsOther = 0;

// prevent name mangling
extern "C" {
    /**
     * Prints the measured FLOPS.
     */
  void printFlops() {
    int rank;
    long long totalLibxsmmFlops;
    long long totalNonZeroFlops;
    long long totalHardwareFlops;
    
    long long nonZeroFlops = g_SeisSolNonZeroFlopsLocal + g_SeisSolNonZeroFlopsNeighbor + g_SeisSolNonZeroFlopsOther;
    long long hardwareFlops = g_SeisSolHardwareFlopsLocal + g_SeisSolHardwareFlopsNeighbor + g_SeisSolHardwareFlopsOther;
    
#ifdef USE_MPI
    long long maxHardwareFlops;

    MPI_Reduce(&libxsmm_num_total_flops, &totalLibxsmmFlops, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&nonZeroFlops, &totalNonZeroFlops, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&hardwareFlops, &totalHardwareFlops, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&hardwareFlops, &maxHardwareFlops, 1, MPI_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  
    double loadImbalance = maxHardwareFlops - ((double) totalHardwareFlops) / size;
    logInfo(rank) << "Load imbalance:            " << loadImbalance;
    logInfo(rank) << "Relative load imbalance:   " << loadImbalance / maxHardwareFlops * 100.0 << "%";
#else
    rank = 0;
    totalLibxsmmFlops = libxsmm_num_total_flops;
    totalNonZeroFlops = nonZeroFlops;
    totalHardwareFlops = hardwareFlops;
#endif

    logInfo(rank) << "Total   measured HW-GFLOP: " << ((double)totalLibxsmmFlops)/1e9;
    logInfo(rank) << "Total calculated HW-GFLOP: " << ((double)totalHardwareFlops)/1e9;
    logInfo(rank) << "Total calculated NZ-GFLOP: " << ((double)totalNonZeroFlops)/1e9;
  }
}
