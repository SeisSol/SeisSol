// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MONITORING_ENERGY_DISCOVERY_H_
#define SEISSOL_SRC_MONITORING_ENERGY_DISCOVERY_H_

#include "Monitoring/Energy/Sensor.h"

#include <memory>
#include <string>
#include <vector>

namespace seissol::monitoring::energy {

/**
 * @brief The sensors of the node, with at most one source for each part, below @p root (/ on a
 * machine).
 *
 * The node, its CPUs, its memory and its accelerators as the HPE Cray power management counters
 * count them; the CPUs and their memory through RAPL where those counters are not there. The GPUs
 * of the node through the management library of their vendor, as found on the PCI bus: NVML for
 * NVIDIA, AMD SMI (else ROCm SMI) for AMD, Level Zero (else the hwmon counters of their driver) for
 * Intel; all of them, whether the run uses them or not. With @p devices false, no GPU library is
 * loaded.
 */
std::vector<std::unique_ptr<Sensor>> nodeSensors(const std::string& root, bool devices = true);

} // namespace seissol::monitoring::energy

#endif // SEISSOL_SRC_MONITORING_ENERGY_DISCOVERY_H_
