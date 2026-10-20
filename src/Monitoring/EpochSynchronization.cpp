// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "EpochSynchronization.h"

#include "Modules/Modules.h"

#include <utils/logger.h>

namespace seissol::monitoring {

void EpochSynchronization::init(double interval) {
  if (interval <= 0) {
    return;
  }
  logInfo() << "Ending a simulation epoch at least every" << interval << "s.";
  setSyncInterval(interval);
  Modules::registerHook(*this, ModuleHook::SynchronizationPoint);
}

} // namespace seissol::monitoring
