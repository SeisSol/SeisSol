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
#include "Monitoring/Counters.h"
#include "Monitoring/Metric.h"

#include <fstream>
#include <map>
#include <string>
#include <vector>

// Floating point operations performed in the matrix kernels.
// Remark: These variables are updated by the matrix kernels (subroutine.cpp) only in debug builds.
// The GEMM generators write these names into the kernels.

// NOLINTNEXTLINE(readability-identifier-naming,cppcoreguidelines-avoid-non-const-global-variables)
extern long long libxsmm_num_total_flops;
// NOLINTNEXTLINE(readability-identifier-naming,cppcoreguidelines-avoid-non-const-global-variables)
extern long long pspamm_num_total_flops;

namespace seissol::monitoring {

struct FlopCounter {
  public:
  void init(const std::string& outputFileNamePrefix, const std::string& runId);

  /**
   * @brief Registers the hardware and non-zero FLOP and the kernel bytes, summed over all metrics,
   * as counters of the run.
   */
  void registerCounters(Counters& counters);

  /**
   * @brief Prints the rates since the start and of the epoch that ended at @p wallTime, and writes
   * them per rank to the FLOP table.
   *
   * Collective.
   */
  void printPerformanceUpdate(double wallTime);

  /**
   * @brief Prints the FLOP and bytes of the whole simulation, which took @p wallTime.
   *
   * Collective. Returns what it printed, for the run summary.
   */
  json::Value printPerformanceSummary(double wallTime) const;

  void incrementMetric(std::size_t handle, const PerformanceEstimate& data);

  [[nodiscard]] std::size_t addMetric(const std::string& name, const std::string& category);

  private:
  std::ofstream out_;
  std::string runId_;

  struct Metric {
    PerformanceEstimate estimate;
    std::string name;
    std::string category;
  };

  std::vector<Metric> metrics_;
  std::map<std::string, std::vector<std::size_t>> metricCategories_;

  PerformanceEstimate previousEstimate_;
  double previousWallTime_ = 0;
};
} // namespace seissol::monitoring

#endif // SEISSOL_SRC_MONITORING_FLOPCOUNTER_H_
