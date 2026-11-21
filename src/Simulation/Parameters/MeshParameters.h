// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SIMULATION_PARAMETERS_MESHPARAMETERS_H_
#define SEISSOL_SRC_SIMULATION_PARAMETERS_MESHPARAMETERS_H_

#include "Core/Geometry/MeshParameters.h"
#include "ParameterReader.h"

namespace seissol::initializer::parameters {

MeshParameters readMeshParameters(ParameterReader* baseReader);
} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_SIMULATION_PARAMETERS_MESHPARAMETERS_H_
