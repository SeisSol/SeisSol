// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SIMULATION_PARAMETERS_MODELPARAMETERS_H_
#define SEISSOL_SRC_SIMULATION_PARAMETERS_MODELPARAMETERS_H_

#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/RunConfigs.h"
#include "Features/Plasticity/PlasticityParameters.h"
#include "Features/TimeMirror/ITMParameters.h"
#include "ParameterReader.h"
#include "WavePropagation/WavePropagationParameters.h"

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

namespace seissol::initializer::parameters {

/// The files of the models of a run, which the composition opens (section `equations`).
struct ModelParameters {
  bool hasBoundaryFile{false};
  // one per simulation, empty where there is none
  std::vector<std::string> boundaryFileNames;
  // one per simulation; the simulations fused into a cell name the same one
  std::vector<std::string> materialFileNames;
  // one per simulation
  std::vector<std::string> plasticityFileNames;
};

/// The configuration the cells of the run compute in: the one named by `configuration` in the
/// section `equations`, or else the one the environment variable SEISSOL_CONFIGURATION names, or
/// else the first one built into the executable (see `defaultConfig`).
ConfigId readConfig(ParameterReader* baseReader);
/// The mesh groups whose cells compute in a configuration other than `config`, by `configmap` in
/// the section `equations`: e.g. "1,2:name;3:other" puts the groups 1 and 2 into the configuration
/// `name`, and the group 3 into `other`.
std::unordered_map<int, ConfigId> readGroupConfigs(ParameterReader* baseReader, ConfigId config);
/// The number of simulations of the run, by `simcount` in the section `equations`; the number of
/// fused simulations of `config` if not given.
std::size_t readSimulationCount(ParameterReader* baseReader, ConfigId config);
/// The configurations of the run and its number of simulations (see the functions above).
RunConfigs readRunConfigs(ParameterReader* baseReader);
ModelParameters readModelParameters(ParameterReader* baseReader, const RunConfigs& runConfigs);
WavePropagationParameters readWavePropagationParameters(ParameterReader* baseReader,
                                                        const RunConfigs& runConfigs);
PlasticityParameters readPlasticityParameters(ParameterReader* baseReader);
ITMParameters readITMParameters(ParameterReader* baseReader);
} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_SIMULATION_PARAMETERS_MODELPARAMETERS_H_
