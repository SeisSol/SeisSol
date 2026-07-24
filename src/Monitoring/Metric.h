// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_MONITORING_METRIC_H_
#define SEISSOL_SRC_MONITORING_METRIC_H_

#include <cstdint>
namespace seissol {

struct PerformanceEstimate {
  std::uint64_t hardwareFlop{0};
  std::uint64_t nonzeroFlop{0};
  std::uint64_t bytes{0};
  std::uint64_t kernelBytes{0};

  constexpr auto operator+(const PerformanceEstimate& other) const -> PerformanceEstimate {
    return PerformanceEstimate{hardwareFlop + other.hardwareFlop,
                               nonzeroFlop + other.nonzeroFlop,
                               bytes + other.bytes,
                               kernelBytes + other.kernelBytes};
  }
};

} // namespace seissol
#endif // SEISSOL_SRC_MONITORING_METRIC_H_
