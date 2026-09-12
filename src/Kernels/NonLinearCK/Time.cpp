// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Time.h"

#include "Initializer/Typedefs.h"

namespace seissol::kernels::solver::nonlinearck {

void Spacetime::setGlobalData(const CompoundGlobalData& global) {
  derivative_.bindGlobals(*global.onHost);
  step_.bindGlobals(*global.onHost);

#ifdef ACL_DEVICE
  deviceDerivative_.bindGlobals(*global.onDevice);
  deviceStep_.bindGlobals(*global.onDevice);
#endif
}

void Time::setGlobalData(const CompoundGlobalData& global) {}

} // namespace seissol::kernels::solver::nonlinearck
