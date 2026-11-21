// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SIMULATION_PARAMETERS_INITIALIZATIONPARAMETERS_H_
#define SEISSOL_SRC_SIMULATION_PARAMETERS_INITIALIZATIONPARAMETERS_H_

#include "Core/Common/ConfigRegistry.h"
#include "Features/InitialConditions/InitializationParameters.h"
#include "ParameterReader.h"

namespace seissol::initializer::parameters {

InitializationParameters readInitializationParameters(ParameterReader* baseReader, ConfigId config);
} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_SIMULATION_PARAMETERS_INITIALIZATIONPARAMETERS_H_
