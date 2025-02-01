// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Rettenberger

#include "Unit.h"
#include <cassert>
#include <cstddef>
#include <fstream>
#include <mpi.h>
#include <ostream>
#include <string>

#include "Numerical/Statistics.h"

// NOLINTNEXTLINE
long long libxsmm_num_total_flops = 0;
// NOLINTNEXTLINE
long long pspamm_num_total_flops = 0;

#include "Parallel/MPI.h"

#include "FlopCounter.h"

#include <utils/logger.h>

namespace seissol::monitoring {

void FlopCounter::init(const std::string& outputFileNamePrefix) {
  const std::string outputFileName = outputFileNamePrefix + "-flops.csv";
  const int rank = seissol::MPI::mpi.rank();
  const auto worldSize = static_cast<size_t>(seissol::MPI::mpi.size());
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
  const auto worldSize = static_cast<size_t>(seissol::MPI::mpi.size());

  const long long newTotalFlops = hardwareFlopsLocal + hardwareFlopsNeighbor + hardwareFlopsOther +
                                  hardwareFlopsDynamicRupture + hardwareFlopsPlasticity;
  const long long diffFlops = newTotalFlops - previousTotalFlops;
  previousTotalFlops = newTotalFlops;

  const double diffTime = wallTime - previousWallTime;
  previousWallTime = wallTime;

  const double accumulatedGflopsPerSecond = newTotalFlops * 1.e-9 / wallTime;
  const double previousGflopsPerSecond = diffFlops * 1.e-9 / diffTime;

  const auto accumulatedGflopsPerSecondOnRanks =
      seissol::MPI::mpi.collect(accumulatedGflopsPerSecond);
  const auto previousGflopsPerSecondOnRanks = seissol::MPI::mpi.collect(previousGflopsPerSecond);

  const auto accumulatedGflopsPerSecondSummary =
      seissol::statistics::Summary(accumulatedGflopsPerSecondOnRanks);
  const auto previousGflopsPerSecondSummary =
      seissol::statistics::Summary(previousGflopsPerSecondOnRanks);

  if (rank == 0) {
    // for now, we calculate everything in GFLOP/s, and switch back to FLOP/s for output only
    logInfo() << "Performance since the start:"
              << UnitFlopPerS.formatPrefix(accumulatedGflopsPerSecondSummary.sum * 1e9).c_str()
              << "(per rank:"
              << UnitFlopPerS
                     .formatPrefix(accumulatedGflopsPerSecondSummary.mean * 1e9,
                                   accumulatedGflopsPerSecondSummary.std * 1e9)
                     .c_str()
              << ")";
    logInfo() << "Performance since last sync point:"
              << UnitFlopPerS.formatPrefix(previousGflopsPerSecondSummary.sum * 1e9).c_str()
              << "(per rank:"
              << UnitFlopPerS
                     .formatPrefix(previousGflopsPerSecondSummary.mean * 1e9,
                                   previousGflopsPerSecondSummary.std * 1e9)
                     .c_str()
              << ")";

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
void FlopCounter::printPerformanceSummary(double wallTime) const {
  enum Counter {
    Libxsmm = 0,
    Pspamm,
    WPNonZeroFlops,
    WPHardwareFlops,
    DRNonZeroFlops,
    DRHardwareFlops,
    PLNonZeroFlops,
    PLHardwareFlops,
    NumCounters
  };

  double flops[NumCounters];

  flops[Libxsmm] = libxsmm_num_total_flops;
  flops[Pspamm] = pspamm_num_total_flops;
  flops[WPNonZeroFlops] = nonZeroFlopsLocal + nonZeroFlopsNeighbor + nonZeroFlopsOther;
  flops[WPHardwareFlops] = hardwareFlopsLocal + hardwareFlopsNeighbor + hardwareFlopsOther;
  flops[DRNonZeroFlops] = nonZeroFlopsDynamicRupture;
  flops[DRHardwareFlops] = hardwareFlopsDynamicRupture;
  flops[PLNonZeroFlops] = nonZeroFlopsPlasticity;
  flops[PLHardwareFlops] = hardwareFlopsPlasticity;

#ifdef USE_MPI
  double totalFlops[NumCounters];
  MPI_Reduce(&flops, &totalFlops, NumCounters, MPI_DOUBLE, MPI_SUM, 0, seissol::MPI::mpi.comm());
#else
  double* totalFlops = &flops[0];
#endif

#ifndef NDEBUG
  logInfo() << "Total    libxsmm HW-FLOP: " << UnitFlop.formatPrefix(totalFlops[Libxsmm]).c_str();
  logInfo() << "Total     pspamm HW-FLOP: " << UnitFlop.formatPrefix(totalFlops[Pspamm]).c_str();
#endif
  const auto totalHardwareFlops =
      totalFlops[WPHardwareFlops] + totalFlops[DRHardwareFlops] + totalFlops[PLHardwareFlops];
  const auto totalNonZeroFlops =
      totalFlops[WPNonZeroFlops] + totalFlops[DRNonZeroFlops] + totalFlops[PLNonZeroFlops];

  const auto percentageUsefulFlops = totalNonZeroFlops / totalHardwareFlops * 100;

  logInfo() << "Total calculated HW-FLOP: " << UnitFlop.formatPrefix(totalHardwareFlops).c_str();
  logInfo() << "Total calculated NZ-FLOP: " << UnitFlop.formatPrefix(totalNonZeroFlops).c_str();
  logInfo() << "NZ part of HW-FLOP:" << percentageUsefulFlops << "%";
  logInfo() << "Total calculated HW-FLOP/s: "
            << UnitFlopPerS
                   .formatPrefix((totalFlops[WPHardwareFlops] + totalFlops[DRHardwareFlops] +
                                  totalFlops[PLHardwareFlops]) /
                                 wallTime)
                   .c_str();
  logInfo() << "Total calculated NZ-FLOP/s: "
            << UnitFlopPerS
                   .formatPrefix((totalFlops[WPNonZeroFlops] + totalFlops[DRNonZeroFlops] +
                                  totalFlops[PLNonZeroFlops]) /
                                 wallTime)
                   .c_str();
  logInfo() << "WP calculated HW-FLOP: "
            << UnitFlop.formatPrefix(totalFlops[WPHardwareFlops]).c_str();
  logInfo() << "WP calculated NZ-FLOP: "
            << UnitFlop.formatPrefix(totalFlops[WPNonZeroFlops]).c_str();
  logInfo() << "DR calculated HW-FLOP: "
            << UnitFlop.formatPrefix(totalFlops[DRHardwareFlops]).c_str();
  logInfo() << "DR calculated NZ-FLOP: "
            << UnitFlop.formatPrefix(totalFlops[DRNonZeroFlops]).c_str();
  logInfo() << "PL calculated HW-FLOP: "
            << UnitFlop.formatPrefix(totalFlops[PLHardwareFlops]).c_str();
  logInfo() << "PL calculated NZ-FLOP: "
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
