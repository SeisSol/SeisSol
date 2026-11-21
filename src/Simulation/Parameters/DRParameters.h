// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SIMULATION_PARAMETERS_DRPARAMETERS_H_
#define SEISSOL_SRC_SIMULATION_PARAMETERS_DRPARAMETERS_H_

#include "Core/Common/ConfigRegistry.h"
#include "Features/DynamicRupture/DRParameters.h"
#include "ParameterReader.h"

namespace seissol::initializer::parameters {

DRParameters readDRParameters(ParameterReader* baseReader, ConfigId config);

} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_SIMULATION_PARAMETERS_DRPARAMETERS_H_
