// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Sensor.h"

#include <functional>
#include <string_view>
#include <vector>

namespace seissol::monitoring::energy {

std::vector<Part> basisParts(const std::function<bool(Part)>& measured) {
  if (measured(Part::Node)) {
    return {Part::Node};
  }
  std::vector<Part> parts;
  for (const auto part : {Part::Cpu, Part::Memory}) {
    if (measured(part)) {
      parts.push_back(part);
    }
  }
  if (measured(Part::Devices)) {
    parts.push_back(Part::Devices);
  } else if (measured(Part::Accelerators)) {
    parts.push_back(Part::Accelerators);
  }
  return parts;
}

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
