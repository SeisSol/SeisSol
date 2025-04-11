// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_INTEGRATIONDATA_H_
#define SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_INTEGRATIONDATA_H_

// #include "Datastructures.h"
#include "Kernels/Precision.h"
#include "generated_code/tensor.h"
#include <Geometry/MeshDefinition.h>

namespace seissol::model {

struct DamageLocalData {
    Vertex localVertices[4];
};
struct DamageNeighborData {};

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_INTEGRATIONDATA_H_
