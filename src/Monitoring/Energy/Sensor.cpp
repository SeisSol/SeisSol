// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Sensor.h"

#include <string_view>

namespace seissol::monitoring::energy {

std::string_view partName(Part part) {
  switch (part) {
  case Part::Node:
    return "node";
  case Part::Cpu:
    return "cpu";
  case Part::Memory:
    return "memory";
  case Part::Accelerators:
    return "accelerators";
  case Part::Devices:
    return "devices";
  }
  return "unknown";
}

double Accumulator::update(double raw) {
  if (!first_.has_value()) {
    first_ = raw;
    last_ = raw;
    return 0;
  }
  if (raw < last_) {
    // a counter with a range ran over it; one without was reset to zero
    offset_ += range_.value_or(last_);
  }
  last_ = raw;
  return raw + offset_ - first_.value();
}

} // namespace seissol::monitoring::energy
