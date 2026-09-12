// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Local.h"

#include "Initializer/Typedefs.h"

namespace seissol::kernels::solver::nonlinearck {

void Local::setGlobalData(const CompoundGlobalData& global) {
  convertToNodal_.bindGlobals(global.onHost);
  convertToModal_.bindGlobals(global.onHost);
  flux_.bindGlobals(global.onHost);
  projectToFace_.bindGlobals(global.onHost);
  rusanov_.bindGlobals(global.onHost);
  faceIntegral_.bindGlobals(global.onHost);
  projectKrnlPrototype_.bindGlobals(global.onHost);
  projectRotatedKrnlPrototype_.bindGlobals(global.onHost);

#ifdef ACL_DEVICE
  deviceConvertToNodal_.bindGlobals(global.onDevice);
  deviceConvertToModal_.bindGlobals(global.onDevice);
  deviceFlux_.bindGlobals(global.onDevice);
  deviceProjectToFace_.bindGlobals(global.onDevice);
  deviceRusanov_.bindGlobals(global.onDevice);
  deviceFaceIntegral_.bindGlobals(global.onDevice);
  deviceProjectRotatedKrnlPrototype_.bindGlobals(global.onDevice);
#endif
}

} // namespace seissol::kernels::solver::nonlinearck
