// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer

#ifndef SEISSOL_SRC_MONITORING_FLOPCOUNTER_H_
#define SEISSOL_SRC_MONITORING_FLOPCOUNTER_H_

#include "Common/Json.h"
#include "Monitoring/Metric.h"

#include <fstream>
#include <map>
#include <mutex>
#include <string>
#include <vector>

// Floating point operations performed in the matrix kernels.
// Remark: These variables are updated by the matrix kernels (subroutine.cpp) only in debug builds.

// NOLINTNEXTLINE
extern long long libxsmm_num_total_flops;
// NOLINTNEXTLINE
extern long long pspamm_num_total_flops;

namespace seissol::monitoring {
//! @brief The rates of one epoch, summed over all ranks.
struct EpochRates {
  double hardwareFlopPerSecond{0};
  double nonzeroFlopPerSecond{0};
  double kernelBytesPerSecond{0};
};

struct FlopCounter {
  public:
  void init(const std::string& outputFileNamePrefix, const std::string& runId);

  /**
   * @brief Prints the rates since the start and of the epoch that ended at @p wallTime.
   *
   * Collective. Returns the rates of the epoch (complete on rank 0 only).
   */
  EpochRates printPerformanceUpdate(double wallTime);

  /**
   * @brief Prints the FLOP and bytes of the whole simulation, which took @p wallTime.
   *
   * Collective. Returns what it printed, for the run summary.
   */
  json::Value printPerformanceSummary(double wallTime) const;

  /**
   * @brief Counts @p data for the metric; from any thread (e.g. that of host work for the device).
   */
  void incrementMetric(std::size_t handle, const PerformanceEstimate& data);

  [[nodiscard]] std::size_t addMetric(const std::string& name, const std::string& category);

  /**
   * @brief What the metrics with one of the @p names have counted on this rank so far.
   */
  [[nodiscard]] PerformanceEstimate estimate(const std::vector<std::string>& names) const;

  private:
  std::ofstream out_;
  std::string runId_;

  struct Metric {
    PerformanceEstimate estimate;
    std::string name;
    std::string category;
  };

  std::vector<Metric> metrics_;
  // for the estimates of the metrics, which work on other threads may increment
  mutable std::mutex mutex_;
  std::map<std::string, std::vector<std::size_t>> metricCategories_;

  PerformanceEstimate previousEstimate_;
  double previousWallTime_ = 0;
};
} // namespace seissol::monitoring

#endif // SEISSOL_SRC_MONITORING_FLOPCOUNTER_H_
