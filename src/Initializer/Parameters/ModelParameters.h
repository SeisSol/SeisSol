// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERS_MODELPARAMETERS_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERS_MODELPARAMETERS_H_

#include "Common/ConfigRegistry.h"
#include "ParameterReader.h"

#include <string>
#include <unordered_map>
#include <vector>

namespace seissol::initializer::parameters {

enum class ReflectionType { BothWaves = 1, BothWavesVelocity, Pwave, Swave };

struct ITMParameters {
  bool itmEnabled{false};
  double itmStartingTime{0.0};
  double itmDuration{0.0};
  double itmVelocityScalingFactor{1.0};
  ReflectionType itmReflectionType{ReflectionType::BothWaves};
  // a script giving the material after a mirror, in place of the scaling of the reflection type
  std::string itmMaterialScript;
};

enum class NumericalFlux { Godunov, Rusanov };

std::string fluxToString(NumericalFlux flux);

struct ModelParameters {
  bool hasBoundaryFile{false};
  bool plasticity{false};
  bool plasticityPointwise{true};
  std::unordered_set<int> plasticityDisabledGroups;
  bool useCellHomogenizedMaterial{true};
  double freqCentral{};
  double freqRatio{1.0};
  double gravitationalAcceleration{};
  double tv{};
  std::string boundaryFileName;
  // the script of the nonlinear Dirichlet boundary (FaceType::NonlinearDirichlet); empty for none
  std::string nonlinearDirichletFileName;
  std::string materialFileName;
  std::vector<std::string> plasticityFileNames;
  ITMParameters itmParameters;
  NumericalFlux flux{NumericalFlux::Godunov};
  NumericalFlux fluxNearFault{NumericalFlux::Godunov};
  // the configuration the cells of the run compute in, unless their mesh group has its own
  ConfigId config{};
  // the mesh groups whose cells compute in a configuration of their own
  std::unordered_map<int, ConfigId> groupConfigs;

  /// The configuration the cells of a mesh group compute in.
  [[nodiscard]] ConfigId configOfGroup(int group) const;
  /// Every configuration cells of the run may compute in: `config` first, then those of the mesh
  /// groups in the order of their ids, each once.
  [[nodiscard]] std::vector<ConfigId> configs() const;
};

/// The configuration the cells of the run compute in: the one named by `configuration` in the
/// section `equations`, or the first one built into the executable.
ConfigId readConfig(ParameterReader* baseReader);
/// The mesh groups whose cells compute in a configuration other than `config`, by `configmap` in
/// the section `equations`: e.g. "1,2:name;3:other" puts the groups 1 and 2 into the configuration
/// `name`, and the group 3 into `other`.
std::unordered_map<int, ConfigId> readGroupConfigs(ParameterReader* baseReader, ConfigId config);
ModelParameters readModelParameters(ParameterReader* baseReader,
                                    ConfigId config,
                                    std::unordered_map<int, ConfigId> groupConfigs);
ITMParameters readITMParameters(ParameterReader* baseReader);
} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERS_MODELPARAMETERS_H_
