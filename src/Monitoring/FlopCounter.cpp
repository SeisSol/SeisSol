/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de,
 *http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
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

#include "Unit.hpp"
#include <cassert>
#include <fstream>

long long libxsmm_num_total_flops = 0;
long long pspamm_num_total_flops = 0;

#include "Parallel/MPI.h"

#include "FlopCounter.hpp"

#include <utils/logger.h>

namespace seissol::monitoring {

void FlopCounter::init(std::string outputFileNamePrefix) {
  const std::string outputFileName = outputFileNamePrefix + "-flops.csv";
  const int rank = seissol::MPI::mpi.rank();
  const int worldSize = seissol::MPI::mpi.size();
  if (rank == 0) {
    out.open(outputFileName);
    out << "time,";
    for (size_t i = 0; i < worldSize - 1; ++i) {
      out << "rank_" << i << "_accumulated,";
      out << "rank_" << i << "_current,";
    }
    out << "rank_" << worldSize - 1 << "_accumulated,";
    out << "rank_" << worldSize - 1 << "_current" << std::endl;
  }
}

void FlopCounter::printPerformanceUpdate(double wallTime) {
  const int rank = seissol::MPI::mpi.rank();
  const int worldSize = seissol::MPI::mpi.size();

  const long long newTotalFlops = hardwareFlopsLocal + hardwareFlopsNeighbor + hardwareFlopsOther +
                                  hardwareFlopsDynamicRupture + hardwareFlopsPlasticity;
  const long long diffFlops = newTotalFlops - previousTotalFlops;
  previousTotalFlops = newTotalFlops;

  const double diffTime = wallTime - previousWallTime;
  previousWallTime = wallTime;

  const double accumulatedGflopsPerSecond = newTotalFlops * 1.e-9 / wallTime;
  const double previousGflopsPerSecond = diffFlops * 1.e-9 / diffTime;

  auto accumulatedGflopsPerSecondOnRanks = seissol::MPI::mpi.collect(accumulatedGflopsPerSecond);
  auto previousGflopsPerSecondOnRanks = seissol::MPI::mpi.collect(previousGflopsPerSecond);

  if (rank == 0) {
    double accumulatedGflopsSum = 0;
    double previousGflopsSum = 0;
#ifdef _OPENMP
#pragma omp simd reduction(+ : accumulatedGflopsSum, previousGflopsSum)
#endif
    for (size_t i = 0; i < worldSize; i++) {
      accumulatedGflopsSum += accumulatedGflopsPerSecondOnRanks[i];
      previousGflopsSum += previousGflopsPerSecondOnRanks[i];
    }
    const auto accumulatedGflopsPerRank = accumulatedGflopsSum / seissol::MPI::mpi.size();
    const auto previousGflopsPerRank = previousGflopsSum / seissol::MPI::mpi.size();

    // for now, we calculate everything in GFLOP/s, and switch back to FLOP/s for output only
    logInfo(rank) << "Performance since the start:"
                  << UnitFlopPerS.formatPrefix(accumulatedGflopsSum * 1e9).c_str() << "(rank 0:"
                  << UnitFlopPerS.formatPrefix(accumulatedGflopsPerSecond * 1e9).c_str()
                  << ", average over ranks:"
                  << UnitFlopPerS.formatPrefix(accumulatedGflopsPerRank * 1e9).c_str() << ")";
    logInfo(rank) << "Performance since last sync point:"
                  << UnitFlopPerS.formatPrefix(previousGflopsSum * 1e9).c_str()
                  << "(rank 0:" << UnitFlopPerS.formatPrefix(previousGflopsPerSecond * 1e9).c_str()
                  << ", average over ranks:"
                  << UnitFlopPerS.formatPrefix(previousGflopsPerRank * 1e9).c_str() << ")";

    out << wallTime << ",";
    for (size_t i = 0; i < worldSize - 1; i++) {
      out << accumulatedGflopsPerSecondOnRanks[i] << ",";
      out << previousGflopsPerSecondOnRanks[i] << ",";
    }
    out << accumulatedGflopsPerSecondOnRanks[worldSize - 1] << ","
        << previousGflopsPerSecondOnRanks[worldSize - 1] << std::endl;
  }
}

/**
 * Prints the measured FLOP/s.
 */
void FlopCounter::printPerformanceSummary(double wallTime) {
  const int rank = seissol::MPI::mpi.rank();

  enum Counter {
    Libxsmm = 0,
    Pspamm,
    WPNonZeroFlops,
    WPHardwareFlops,
    DRNonZeroFlops,
    DRHardwareFlops,
    PLNonZeroFlops,
    PLHardwareFlops,
    NUM_COUNTERS
  };

  double flops[NUM_COUNTERS];

  flops[Libxsmm] = libxsmm_num_total_flops;
  flops[Pspamm] = pspamm_num_total_flops;
  flops[WPNonZeroFlops] = nonZeroFlopsLocal + nonZeroFlopsNeighbor + nonZeroFlopsOther;
  flops[WPHardwareFlops] = hardwareFlopsLocal + hardwareFlopsNeighbor + hardwareFlopsOther;
  flops[DRNonZeroFlops] = nonZeroFlopsDynamicRupture;
  flops[DRHardwareFlops] = hardwareFlopsDynamicRupture;
  flops[PLNonZeroFlops] = nonZeroFlopsPlasticity;
  flops[PLHardwareFlops] = hardwareFlopsPlasticity;

#ifdef USE_MPI
  double totalFlops[NUM_COUNTERS];
  MPI_Reduce(&flops, &totalFlops, NUM_COUNTERS, MPI_DOUBLE, MPI_SUM, 0, seissol::MPI::mpi.comm());
#else
  double* totalFlops = &flops[0];
#endif

#ifndef NDEBUG
  logInfo(rank) << "Total    libxsmm HW-FLOP: "
                << UnitFlop.formatPrefix(totalFlops[Libxsmm]).c_str();
  logInfo(rank) << "Total     pspamm HW-FLOP: "
                << UnitFlop.formatPrefix(totalFlops[Pspamm]).c_str();
#endif
  logInfo(rank) << "Total calculated HW-FLOP: "
                << UnitFlop
                       .formatPrefix(totalFlops[WPHardwareFlops] + totalFlops[DRHardwareFlops] +
                                     totalFlops[PLHardwareFlops])
                       .c_str();
  logInfo(rank) << "Total calculated NZ-FLOP: "
                << UnitFlop
                       .formatPrefix(totalFlops[WPNonZeroFlops] + totalFlops[DRNonZeroFlops] +
                                     totalFlops[PLNonZeroFlops])
                       .c_str();
  logInfo(rank) << "Total calculated HW-FLOP/s: "
                << UnitFlopPerS
                       .formatPrefix((totalFlops[WPHardwareFlops] + totalFlops[DRHardwareFlops] +
                                      totalFlops[PLHardwareFlops]) /
                                     wallTime)
                       .c_str();
  logInfo(rank) << "Total calculated NZ-FLOP/s: "
                << UnitFlopPerS
                       .formatPrefix((totalFlops[WPNonZeroFlops] + totalFlops[DRNonZeroFlops] +
                                      totalFlops[PLNonZeroFlops]) /
                                     wallTime)
                       .c_str();
  logInfo(rank) << "WP calculated HW-FLOP: "
                << UnitFlop.formatPrefix(totalFlops[WPHardwareFlops]).c_str();
  logInfo(rank) << "WP calculated NZ-FLOP: "
                << UnitFlop.formatPrefix(totalFlops[WPNonZeroFlops]).c_str();
  logInfo(rank) << "DR calculated HW-FLOP: "
                << UnitFlop.formatPrefix(totalFlops[DRHardwareFlops]).c_str();
  logInfo(rank) << "DR calculated NZ-FLOP: "
                << UnitFlop.formatPrefix(totalFlops[DRNonZeroFlops]).c_str();
  logInfo(rank) << "PL calculated HW-FLOP: "
                << UnitFlop.formatPrefix(totalFlops[PLHardwareFlops]).c_str();
  logInfo(rank) << "PL calculated NZ-FLOP: "
                << UnitFlop.formatPrefix(totalFlops[PLNonZeroFlops]).c_str();
}
void FlopCounter::incrementNonZeroFlopsLocal(long long update) {
  assert(update >= 0);
  nonZeroFlopsLocal += update;
}
void FlopCounter::incrementHardwareFlopsLocal(long long update) {
  assert(update >= 0);
  hardwareFlopsLocal += update;
}
void FlopCounter::incrementNonZeroFlopsNeighbor(long long update) {
  assert(update >= 0);
  nonZeroFlopsNeighbor += update;
}
void FlopCounter::incrementHardwareFlopsNeighbor(long long update) {
  assert(update >= 0);
  hardwareFlopsNeighbor += update;
}
void FlopCounter::incrementNonZeroFlopsOther(long long update) {
  assert(update >= 0);
  nonZeroFlopsOther += update;
}
void FlopCounter::incrementHardwareFlopsOther(long long update) {
  assert(update >= 0);
  hardwareFlopsOther += update;
}
void FlopCounter::incrementNonZeroFlopsDynamicRupture(long long update) {
  assert(update >= 0);
  nonZeroFlopsDynamicRupture += update;
}
void FlopCounter::incrementHardwareFlopsDynamicRupture(long long update) {
  assert(update >= 0);
  hardwareFlopsDynamicRupture += update;
}
void FlopCounter::incrementNonZeroFlopsPlasticity(long long update) {
  assert(update >= 0);
  nonZeroFlopsPlasticity += update;
}
void FlopCounter::incrementHardwareFlopsPlasticity(long long update) {
  assert(update >= 0);
  hardwareFlopsPlasticity += update;
}
} // namespace seissol::monitoring
