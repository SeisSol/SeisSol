// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SIMULATION_PARAMETERS_SOURCEPARAMETERS_H_
#define SEISSOL_SRC_SIMULATION_PARAMETERS_SOURCEPARAMETERS_H_

#include "Features/SourceTerm/SourceParameters.h"
#include "ParameterReader.h"

namespace seissol::initializer::parameters {

SourceParameters readSourceParameters(ParameterReader* baseReader);
} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_SIMULATION_PARAMETERS_SOURCEPARAMETERS_H_
