// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer

#ifndef SEISSOL_SRC_MONITORING_FLOPCOUNTER_H_
#define SEISSOL_SRC_MONITORING_FLOPCOUNTER_H_

#include <fstream>

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
  void incrementNonZeroFlopsLocal(long long update);
  void incrementHardwareFlopsLocal(long long update);
  void incrementNonZeroFlopsNeighbor(long long update);
  void incrementHardwareFlopsNeighbor(long long update);
  void incrementNonZeroFlopsOther(long long update);
  void incrementHardwareFlopsOther(long long update);
  void incrementNonZeroFlopsDynamicRupture(long long update);
  void incrementHardwareFlopsDynamicRupture(long long update);
  void incrementNonZeroFlopsPlasticity(long long update);
  void incrementHardwareFlopsPlasticity(long long update);

  private:
  std::ofstream out_;
  long long previousTotalHWFlops_ = 0;
  long long previousTotalNZFlops_ = 0;
  double previousWallTime_ = 0;
  // global variables for summing-up SeisSol internal counters
  long long nonZeroFlopsLocal_ = 0;
  long long hardwareFlopsLocal_ = 0;
  long long nonZeroFlopsNeighbor_ = 0;
  long long hardwareFlopsNeighbor_ = 0;
  long long nonZeroFlopsOther_ = 0;
  long long hardwareFlopsOther_ = 0;
  long long nonZeroFlopsDynamicRupture_ = 0;
  long long hardwareFlopsDynamicRupture_ = 0;
  long long nonZeroFlopsPlasticity_ = 0;
  long long hardwareFlopsPlasticity_ = 0;
};
} // namespace seissol::monitoring

#endif // SEISSOL_SRC_MONITORING_FLOPCOUNTER_H_
