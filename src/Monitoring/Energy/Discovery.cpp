// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Discovery.h"

#include "Monitoring/Energy/Sensor.h"
#include "Monitoring/Energy/Sysfs.h"
#include "Monitoring/Energy/Vendors.h"

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace seissol::monitoring::energy {

namespace {
void append(std::vector<std::unique_ptr<Sensor>>& sensors,
            std::vector<std::unique_ptr<Sensor>> more) {
  for (auto& sensor : more) {
    sensors.push_back(std::move(sensor));
  }
}

bool measures(const std::vector<std::unique_ptr<Sensor>>& sensors, Part part) {
  return std::any_of(
      sensors.begin(), sensors.end(), [&](const auto& sensor) { return sensor->part() == part; });
}
} // namespace

std::vector<std::unique_ptr<Sensor>> nodeSensors(const std::string& root, bool devices) {
  auto sensors = crayPmSensors(root);

  // RAPL for what pm_counters do not count, so that no part is counted twice; all packages count
  const bool cpu = measures(sensors, Part::Cpu);
  const bool memory = measures(sensors, Part::Memory);
  for (auto& sensor : powercapSensors(root)) {
    if ((sensor->part() == Part::Cpu && !cpu) || (sensor->part() == Part::Memory && !memory)) {
      sensors.push_back(std::move(sensor));
    }
  }

  if (!devices) {
    return sensors;
  }
  for (const auto vendor : gpuVendors(root)) {
    switch (vendor) {
    case Vendor::Nvidia:
      append(sensors, nvmlSensors());
      break;
    case Vendor::Amd: {
      auto amd = amdSmiSensors();
      append(sensors, amd.empty() ? rocmSmiSensors() : std::move(amd));
      break;
    }
    case Vendor::Intel: {
      auto intel = levelZeroSensors();
      append(sensors, intel.empty() ? hwmonGpuSensors(root) : std::move(intel));
      break;
    }
    }
  }
  return sensors;
}

} // namespace seissol::monitoring::energy
