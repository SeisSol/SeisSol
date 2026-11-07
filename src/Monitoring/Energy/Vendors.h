// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MONITORING_ENERGY_VENDORS_H_
#define SEISSOL_SRC_MONITORING_ENERGY_VENDORS_H_

#include "Monitoring/Energy/Sensor.h"

#include <memory>
#include <vector>

// The energy counters of GPUs, from the management libraries of their vendors, which are loaded at
// runtime. Each function gives a sensor per GPU of the node that its library can read the energy
// of, and none if the library is not there. The declarations of the libraries are those of their
// headers (nvml.h, amd_smi/amdsmi.h, rocm_smi/rocm_smi.h, zes_api.h), for the few functions used.
namespace seissol::monitoring::energy {

//! @brief NVIDIA GPUs, through NVML: the energy since the driver was loaded, in millijoules.
std::vector<std::unique_ptr<Sensor>> nvmlSensors();

//! @brief AMD GPUs, through AMD SMI: an accumulator of a resolution given in microjoules.
std::vector<std::unique_ptr<Sensor>> amdSmiSensors();

//! @brief AMD GPUs, through ROCm SMI, which ROCm installations before AMD SMI have.
std::vector<std::unique_ptr<Sensor>> rocmSmiSensors();

/**
 * @brief Intel GPUs, through the Sysman API of Level Zero (zesInit): the energy of the whole card
 * where the driver gives it, else of the package.
 */
std::vector<std::unique_ptr<Sensor>> levelZeroSensors();

} // namespace seissol::monitoring::energy

#endif // SEISSOL_SRC_MONITORING_ENERGY_VENDORS_H_
