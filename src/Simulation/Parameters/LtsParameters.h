// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SIMULATION_PARAMETERS_LTSPARAMETERS_H_
#define SEISSOL_SRC_SIMULATION_PARAMETERS_LTSPARAMETERS_H_

#include "Core/Common/ConfigRegistry.h"
#include "Core/TimeStepping/LtsParameters.h"
#include "ParameterReader.h"

#include <string>
#include <vector>

namespace seissol::initializer::parameters {

AutoMergeCostBaseline parseAutoMergeCostBaseline(std::string str);

LtsClusteringSearch parseLtsClusteringSearch(std::string str);

LtsParameters readLtsParameters(ParameterReader* baseReader);
/// The time stepping of a run whose cells compute in the configurations `configs`.
TimeSteppingParameters readTimeSteppingParameters(ParameterReader* baseReader,
                                                  const std::vector<ConfigId>& configs);

} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_SIMULATION_PARAMETERS_LTSPARAMETERS_H_
