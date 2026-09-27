// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_INITIALIZER_MODEL_DYNAMICRUPTUREMATRICES_H_
#define SEISSOL_SRC_INITIALIZER_MODEL_DYNAMICRUPTUREMATRICES_H_

#include "Geometry/MeshReader.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"

namespace seissol::initializer {

/// \param global the constant pool, for the kernel that reads a material
/// varying along a face at the quadrature points of the fault
void initializeDynamicRuptureMatrices(const seissol::geometry::MeshReader& meshReader,
                                      LTS::Storage& ltsStorage,
                                      const LTS::Backmap& backmap,
                                      DynamicRupture::Storage& drStorage,
                                      const GlobalData& global);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_MODEL_DYNAMICRUPTUREMATRICES_H_
