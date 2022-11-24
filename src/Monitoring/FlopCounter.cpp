/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2013-2016, SeisSol Group
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

#include <fstream>

#include "Parallel/MPI.h"

#include "FlopCounter.hpp"

#include <utils/logger.h>

namespace seissol::monitoring {

void FlopCounter::init(std::string outputFileNamePrefix) {
  const std::string outputFileName = outputFileNamePrefix + "-flops.csv";
  const int worldSize = seissol::MPI::mpi.size();
  out.open(outputFileName);
  out << "time,";
  for (size_t i = 0; i < worldSize - 1; ++i) {
    out << "rank_" << i << ",";
  }
  out << "rank_" << worldSize - 1 << std::endl;
}

void FlopCounter::printPerformance(double wallTime) {
  const int rank = seissol::MPI::mpi.rank();
  const int worldSize = seissol::MPI::mpi.size();
  const long long flops = g_SeisSolHardwareFlopsLocal
                    + g_SeisSolHardwareFlopsNeighbor
                    + g_SeisSolHardwareFlopsOther
                    + g_SeisSolHardwareFlopsDynamicRupture
                    + g_SeisSolHardwareFlopsPlasticity;
  const double gflopsPerSecond = flops * 1.e-9 / wallTime;

  double gflopsPerSecondOnRanks[worldSize];
  MPI_Gather(&gflopsPerSecond,
             1,
             MPI_DOUBLE,
             gflopsPerSecondOnRanks,
             1,
             MPI_DOUBLE,
             0,
             seissol::MPI::mpi.comm());

  double flopsSum = 0;
  MPI_Reduce(
      &gflopsPerSecond,
      &flopsSum,
      1,
      MPI_DOUBLE,
      MPI_SUM,
      0,
      seissol::MPI::mpi.comm()
  );
  if (rank == 0) {
    const auto flopsPerRank = flopsSum / seissol::MPI::mpi.size();
    logInfo(rank) << flopsSum * 1.e-3  << "TFLOPS"
    << "(rank 0:" << gflopsPerSecond << "GFLOPS, average over ranks:" << flopsPerRank << "GFLOPS)";
    std::ofstream out;
    out.open("flops.csv");
    out << wallTime << ",";
    for (size_t i = 0; i < worldSize - 1; i++) {
      out << gflopsPerSecondOnRanks[i] << ",";
    }
    out << gflopsPerSecondOnRanks[worldSize - 1] << std::endl;
  }
}
  
/**
 * Prints the measured FLOPS.
 */
void FlopCounter::printFlops() {
  const int rank = seissol::MPI::mpi.rank();

  enum Counter {
    Libxsmm = 0,
    WPNonZeroFlops,
    WPHardwareFlops,
    DRNonZeroFlops,
    DRHardwareFlops,
    PLNonZeroFlops,
    PLHardwareFlops,
    NUM_COUNTERS
  };

  double flops[NUM_COUNTERS];

  flops[Libxsmm]          = libxsmm_num_total_flops;
  flops[WPNonZeroFlops]   = g_SeisSolNonZeroFlopsLocal + g_SeisSolNonZeroFlopsNeighbor + g_SeisSolNonZeroFlopsOther;
  flops[WPHardwareFlops]  = g_SeisSolHardwareFlopsLocal + g_SeisSolHardwareFlopsNeighbor + g_SeisSolHardwareFlopsOther;
  flops[DRNonZeroFlops]   = g_SeisSolNonZeroFlopsDynamicRupture;
  flops[DRHardwareFlops]  = g_SeisSolHardwareFlopsDynamicRupture;
  flops[PLNonZeroFlops]   = g_SeisSolNonZeroFlopsPlasticity;
  flops[PLHardwareFlops]  = g_SeisSolHardwareFlopsPlasticity;

#ifdef USE_MPI
  double totalFlops[NUM_COUNTERS];
  MPI_Reduce(&flops, &totalFlops, NUM_COUNTERS, MPI_DOUBLE, MPI_SUM, 0, seissol::MPI::mpi.comm());
#else
  double* totalFlops = &flops[0];
#endif

  logInfo(rank) << "Total   measured HW-GFLOP: " << totalFlops[Libxsmm] * 1.e-9;
  logInfo(rank) << "Total calculated HW-GFLOP: " << (totalFlops[WPHardwareFlops] + totalFlops[DRHardwareFlops] + totalFlops[PLHardwareFlops]) * 1.e-9;
  logInfo(rank) << "Total calculated NZ-GFLOP: " << (totalFlops[WPNonZeroFlops]  + totalFlops[DRNonZeroFlops]  + totalFlops[PLNonZeroFlops] ) * 1.e-9;
  logInfo(rank) << "WP calculated HW-GFLOP: " << (totalFlops[WPHardwareFlops]) * 1.e-9;
  logInfo(rank) << "WP calculated NZ-GFLOP: " << (totalFlops[WPNonZeroFlops])  * 1.e-9;
  logInfo(rank) << "DR calculated HW-GFLOP: " << (totalFlops[DRHardwareFlops]) * 1.e-9;
  logInfo(rank) << "DR calculated NZ-GFLOP: " << (totalFlops[DRNonZeroFlops])  * 1.e-9;
  logInfo(rank) << "PL calculated HW-GFLOP: " << (totalFlops[PLHardwareFlops]) * 1.e-9;
  logInfo(rank) << "PL calculated NZ-GFLOP: " << (totalFlops[PLNonZeroFlops])  * 1.e-9;
}
}
