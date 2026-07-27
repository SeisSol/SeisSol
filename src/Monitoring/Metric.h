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

  constexpr auto operator+=(const PerformanceEstimate& other) -> PerformanceEstimate& {
    *this = *this + other;
    return *this;
  }

  constexpr auto operator*(std::size_t other) const -> PerformanceEstimate {
    return PerformanceEstimate{
        hardwareFlop * other, nonzeroFlop * other, bytes * other, kernelBytes * other};
  }

  constexpr auto operator*=(std::size_t other) -> PerformanceEstimate& {
    *this = *this * other;
    return *this;
  }

  template <typename T, typename... Args>
  static auto fromKernel(Args... kargs) -> PerformanceEstimate {
    return PerformanceEstimate{T::hardwareFlops(kargs...),
                               T::nonZeroFlops(kargs...),
                               0,
                               T::outboundBytes(kargs...) + T::inboundBytes(kargs...)};
  }
};

} // namespace seissol
#endif // SEISSOL_SRC_MONITORING_METRIC_H_
