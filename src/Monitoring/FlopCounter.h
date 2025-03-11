// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
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
  std::ofstream out;
  long long previousTotalHWFlops = 0;
  long long previousTotalNZFlops = 0;
  double previousWallTime = 0;
  // global variables for summing-up SeisSol internal counters
  long long nonZeroFlopsLocal = 0;
  long long hardwareFlopsLocal = 0;
  long long nonZeroFlopsNeighbor = 0;
  long long hardwareFlopsNeighbor = 0;
  long long nonZeroFlopsOther = 0;
  long long hardwareFlopsOther = 0;
  long long nonZeroFlopsDynamicRupture = 0;
  long long hardwareFlopsDynamicRupture = 0;
  long long nonZeroFlopsPlasticity = 0;
  long long hardwareFlopsPlasticity = 0;
};
} // namespace seissol::monitoring

#endif // SEISSOL_SRC_MONITORING_FLOPCOUNTER_H_
