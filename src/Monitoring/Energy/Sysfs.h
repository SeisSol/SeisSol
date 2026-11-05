// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MONITORING_ENERGY_SYSFS_H_
#define SEISSOL_SRC_MONITORING_ENERGY_SYSFS_H_

#include "Monitoring/Energy/Sensor.h"

#include <memory>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <vector>

// The sources of energy that the kernel exposes as files. Each function takes the directory that
// /sys is below (/ on a machine, another directory in the tests).
namespace seissol::monitoring::energy {

//! @brief A value of a file of the HPE Cray power management counters.
struct CrayPmValue {
  double value{0};
  //! J or W
  std::string unit;
  //! the time of the value, in seconds
  std::optional<double> seconds;
};

//! @brief Parses a file of /sys/cray/pm_counters: '<value> <unit> <timestamp> us'.
std::optional<CrayPmValue> parseCrayPm(std::string_view text);

/**
 * @brief The energy counters of an HPE Cray EX node: the node (energy), its CPUs (cpu_energy), its
 * memory (memory_energy) and each of its accelerators (accel<N>_energy).
 */
std::vector<std::unique_ptr<Sensor>> crayPmSensors(const std::string& root);

//! @brief The power cap of an HPE Cray EX node, in watts, if it has one.
std::optional<double> crayPmPowerCap(const std::string& root);

/**
 * @brief The RAPL zones of the CPU packages and of their memory, of those that can be read.
 *
 * On most kernels, only root may read them.
 */
std::vector<std::unique_ptr<Sensor>> powercapSensors(const std::string& root);

//! @brief The vendors of GPUs.
enum class Vendor { Nvidia, Amd, Intel };

std::string_view vendorName(Vendor vendor);

//! @brief The vendors of the GPUs on the PCI bus of the node.
std::set<Vendor> gpuVendors(const std::string& root);

//! @brief The energy counters that the drivers of Intel GPUs (i915, xe) give in hwmon.
std::vector<std::unique_ptr<Sensor>> hwmonGpuSensors(const std::string& root);

//! @brief The content of the file at @p path, without the whitespace around it.
std::optional<std::string> readSysfsFile(const std::string& path);

} // namespace seissol::monitoring::energy

#endif // SEISSOL_SRC_MONITORING_ENERGY_SYSFS_H_
