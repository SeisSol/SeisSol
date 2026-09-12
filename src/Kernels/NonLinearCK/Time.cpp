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
  derivative_.bindGlobals(global.onHost);
  convertToNodal_.bindGlobals(global.onHost);
  invariants_.bindGlobals(global.onHost);
  stress_.bindGlobals(global.onHost);
  flux_.bindGlobals(global.onHost);
  cellState_.bindGlobals(global.onHost);
  source_.bindGlobals(global.onHost);
  projectDerivativeToNodalBoundaryRotated_.bindGlobals(global.onHost);

#ifdef ACL_DEVICE
  deviceDerivative_.bindGlobals(global.onDevice);
  deviceConvertToNodal_.bindGlobals(global.onDevice);
  deviceInvariants_.bindGlobals(global.onDevice);
  deviceStress_.bindGlobals(global.onDevice);
  deviceFlux_.bindGlobals(global.onDevice);
  deviceCellState_.bindGlobals(global.onDevice);
  deviceSource_.bindGlobals(global.onDevice);
  deviceDerivativeToNodalBoundaryRotated_.bindGlobals(global.onDevice);
#endif
}

void Time::setGlobalData(const CompoundGlobalData& global) {}

} // namespace seissol::kernels::solver::nonlinearck
