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
  projectToFace_.bindGlobals(*global.onHost);
  rusanov_.bindGlobals(*global.onHost);
  faceIntegral_.bindGlobals(*global.onHost);

#ifdef ACL_DEVICE
  deviceProjectToFace_.bindGlobals(*global.onDevice);
  deviceRusanov_.bindGlobals(*global.onDevice);
  deviceFaceIntegral_.bindGlobals(*global.onDevice);
#endif
}

} // namespace seissol::kernels::solver::nonlinearck
