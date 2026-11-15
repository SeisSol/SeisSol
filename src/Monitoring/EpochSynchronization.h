// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MONITORING_EPOCHSYNCHRONIZATION_H_
#define SEISSOL_SRC_MONITORING_EPOCHSYNCHRONIZATION_H_

#include "Memory/HostData.h"
#include "Modules/Module.h"
#include "Modules/Modules.h"

namespace seissol::monitoring {

/**
 * @brief Synchronization points at the multiples of an interval of simulated time, which write
 * nothing.
 *
 * Each ends a simulation epoch, which the simulator reports, as it does at every synchronization
 * point: in the log, the run summary and the epoch table. A run whose outputs synchronize rarely
 * or not at all, such as a benchmark, reports its progress and performance as it goes this way.
 */
class EpochSynchronization : public Module {
  public:
  /**
   * @brief Synchronizes at the multiples of @p interval ; not at all if it is not positive.
   *
   * To be called before the simulation starts; registers with @p modules for it.
   */
  void init(Modules& modules, double interval);

  [[nodiscard]] initializer::HostData hostData() const override { return {}; }

  [[nodiscard]] bool processesAsynchronously() const override { return true; }
};

} // namespace seissol::monitoring

#endif // SEISSOL_SRC_MONITORING_EPOCHSYNCHRONIZATION_H_
