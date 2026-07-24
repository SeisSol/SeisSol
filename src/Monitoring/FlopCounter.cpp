// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Rettenberger

#include "FlopCounter.h"

#include "Monitoring/Metric.h"
#include "Numerical/Statistics.h"
#include "Parallel/MPI.h"
#include "Unit.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <fstream>
#include <mpi.h>
#include <ostream>
#include <string>
#include <utils/logger.h>

// NOLINTNEXTLINE
long long libxsmm_num_total_flops = 0;
// NOLINTNEXTLINE
long long pspamm_num_total_flops = 0;

namespace seissol::monitoring {

void FlopCounter::init(const std::string& outputFileNamePrefix) {
  const std::string outputFileName = outputFileNamePrefix + "-flops.csv";
  const int rank = seissol::Mpi::mpi.rank();
  const auto worldSize = static_cast<size_t>(seissol::Mpi::mpi.size());
  if (rank == 0) {
    out_.open(outputFileName);
    out_ << "time,";
    const auto datasetHeaders = [&](const std::string& suffix) {
      for (size_t i = 0; i < worldSize; ++i) {
        out_ << "rank_" << i << "_" << suffix << ",";
      }
    };
    datasetHeaders("hw_accumulated");
    datasetHeaders("hw_epoch");
    datasetHeaders("nz_accumulated");
    datasetHeaders("nz_epoch");
    out_ << std::endl;
  }
}

void FlopCounter::printPerformanceUpdate(double wallTime) {
  const int rank = seissol::Mpi::mpi.rank();
  const auto worldSize = static_cast<size_t>(seissol::Mpi::mpi.size());

  auto estimate = PerformanceEstimate{};
  for (const auto& metric : metrics_) {
    estimate = estimate + metric.estimate;
  }

  const auto diffHWFlops = estimate.hardwareFlop - previousEstimate_.hardwareFlop;
  const auto diffNZFlops = estimate.nonzeroFlop - previousEstimate_.nonzeroFlop;

  previousEstimate_ = estimate;

  const double diffTime = wallTime - previousWallTime_;
  previousWallTime_ = wallTime;

  const double accumulatedHWGflopsPerSecond = estimate.hardwareFlop * 1.e-9 / wallTime;
  const double accumulatedNZGflopsPerSecond = estimate.nonzeroFlop * 1.e-9 / wallTime;
  const double previousHWGflopsPerSecond = diffHWFlops * 1.e-9 / diffTime;
  const double previousNZGflopsPerSecond = diffNZFlops * 1.e-9 / diffTime;

  if (rank == 0) {
    out_ << wallTime << ",";
  }

  const auto handleFlopsDataset = [&](auto local, const std::string& message) {
    const auto localOnRanks = seissol::Mpi::mpi.collect(local);
    const auto localSummary = seissol::statistics::Summary(localOnRanks);

    if (rank == 0) {
      // for now, we calculate everything in GFLOP/s, and switch back to FLOP/s for output only
      logInfo()
          << message.c_str() << UnitFlopPerS.formatPrefix(localSummary.sum * 1e9).c_str()
          << "(per rank:"
          << UnitFlopPerS.formatPrefix(localSummary.mean * 1e9, localSummary.std * 1e9).c_str()
          << ")";
      for (size_t i = 0; i < worldSize; i++) {
        out_ << localOnRanks[i] << ",";
      }
    }
  };

  // make sure to keep it the same width. Otherwise, if may become confusing
  handleFlopsDataset(accumulatedHWGflopsPerSecond, "HW-FLOP/s since start:");
  handleFlopsDataset(previousHWGflopsPerSecond, "HW-FLOP/s last epoch: ");
  handleFlopsDataset(accumulatedNZGflopsPerSecond, "NZ-FLOP/s since start:");
  handleFlopsDataset(previousNZGflopsPerSecond, "NZ-FLOP/s last epoch: ");

  out_ << std::endl;
}

/**
 * Prints the measured FLOP/s.
 */
void FlopCounter::printPerformanceSummary(double wallTime) const {
  // LIBXSMM + PSpaMM
  constexpr std::size_t CodegenFlops = 2;
  constexpr std::size_t BroadcastEntries = 2;

  std::vector<double> flops{};

  flops.push_back(libxsmm_num_total_flops);
  flops.push_back(pspamm_num_total_flops);

  // total first; then by category
  {
    auto estimate = PerformanceEstimate{};
    for (const auto& metric : metrics_) {
      estimate = estimate + metric.estimate;
    }
    flops.push_back(estimate.nonzeroFlop);
    flops.push_back(estimate.hardwareFlop);
  }

  for (const auto& [_, handles] : metricCategories_) {
    auto estimate = PerformanceEstimate{};
    for (const auto& handle : handles) {
      estimate = estimate + metrics_[handle].estimate;
    }
    flops.push_back(estimate.nonzeroFlop);
    flops.push_back(estimate.hardwareFlop);
  }

  MPI_Allreduce(
      MPI_IN_PLACE, flops.data(), flops.size(), MPI_DOUBLE, MPI_SUM, seissol::Mpi::mpi.comm());

#ifndef NDEBUG
  logInfo() << "Total    libxsmm HW-FLOP: " << UnitFlop.formatPrefix(flops[0]).c_str();
  logInfo() << "Total     pspamm HW-FLOP: " << UnitFlop.formatPrefix(flops[1]).c_str();
#endif
  const auto totalNonZeroFlops = flops[CodegenFlops + 0];
  const auto totalHardwareFlops = flops[CodegenFlops + 1];

  const auto percentageUsefulFlops = totalNonZeroFlops / totalHardwareFlops * 100;

  logInfo() << "Total calculated HW-FLOP: " << UnitFlop.formatPrefix(totalHardwareFlops).c_str();
  logInfo() << "Total calculated NZ-FLOP: " << UnitFlop.formatPrefix(totalNonZeroFlops).c_str();
  logInfo() << "NZ part of HW-FLOP:" << percentageUsefulFlops << "%";
  logInfo() << "Total calculated HW-FLOP/s: "
            << UnitFlopPerS.formatPrefix((totalHardwareFlops) / wallTime).c_str();
  logInfo() << "Total calculated NZ-FLOP/s: "
            << UnitFlopPerS.formatPrefix((totalNonZeroFlops) / wallTime).c_str();

  std::size_t i = 0;
  for (const auto& [catname, _] : metricCategories_) {
    logInfo()
        << catname.c_str() << "calculated HW-FLOP: "
        << UnitFlop.formatPrefix(flops[CodegenFlops + (i + 1) * BroadcastEntries + 1]).c_str();
    logInfo()
        << catname.c_str() << "calculated NZ-FLOP: "
        << UnitFlop.formatPrefix(flops[CodegenFlops + (i + 1) * BroadcastEntries + 0]).c_str();
    ++i;
  }
}

void FlopCounter::incrementMetric(std::size_t handle, const PerformanceEstimate& data) {
  metrics_[handle].estimate = metrics_[handle].estimate + data;
}

std::size_t FlopCounter::addMetric(const std::string& name, const std::string& category) {
  const auto handle = metrics_.size();

  metrics_.emplace_back(Metric{PerformanceEstimate{}, name, category});

  metricCategories_[category].push_back(handle);

  return handle;
}

} // namespace seissol::monitoring
