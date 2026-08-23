// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer

#ifndef SEISSOL_SRC_MONITORING_FLOPCOUNTER_H_
#define SEISSOL_SRC_MONITORING_FLOPCOUNTER_H_

#include "Monitoring/Metric.h"

#include <fstream>
#include <map>
#include <vector>

// Floating point operations performed in the matrix kernels.
// Remark: These variables are updated by the matrix kernels (subroutine.cpp) only in debug builds.

// NOLINTNEXTLINE
extern long long libxsmm_num_total_flops;
// NOLINTNEXTLINE
extern long long pspamm_num_total_flops;

namespace seissol::monitoring {
struct FlopCounter {
  public:
  void init(const std::string& outputFileNamePrefix);
  void printPerformanceUpdate(double wallTime);
  void printPerformanceSummary(double wallTime) const;

  void incrementMetric(std::size_t handle, const PerformanceEstimate& data);

  [[nodiscard]] std::size_t addMetric(const std::string& name, const std::string& category);

  private:
  std::ofstream out_;

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
