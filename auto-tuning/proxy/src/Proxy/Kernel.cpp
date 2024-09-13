// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include "Kernel.h"

namespace seissol::proxy {

auto PerformanceEstimate::operator+(const PerformanceEstimate& other) -> PerformanceEstimate {
  return PerformanceEstimate{
      hardwareFlop + other.hardwareFlop, nonzeroFlop + other.nonzeroFlop, bytes + other.bytes};
}

} // namespace seissol::proxy
