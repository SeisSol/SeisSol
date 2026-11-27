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

#include <array>
#include <cstddef>
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
};

enum class NumericalFlux { Godunov, Rusanov, Centered, Mixed };

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
  // one per simulation, empty where there is none
  std::vector<std::string> boundaryFileNames;
  // one per simulation; the simulations fused into a cell name the same one
  std::vector<std::string> materialFileNames;
  // one per simulation
  std::vector<std::string> plasticityFileNames;
  ITMParameters itmParameters;
  NumericalFlux flux{NumericalFlux::Godunov};
  NumericalFlux fluxNearFault{NumericalFlux::Godunov};
  // the flux "mixed": the share of the Godunov flux (the rest is the centered flux) and the
  // factor of the dissipation of the Rusanov flux, away from the fault and near it
  double fluxRPart{0};
  double fluxNearFaultRPart{0};
  double fluxGPart{0};
  double fluxNearFaultGPart{0};
  // the shares of the upwind damping in the Godunov part of the flux (see
  // model::scaleRiemannDamping): of the velocity jumps (normal, tangential) and of the traction
  // jumps (normal, tangential); all 1 is the Godunov flux
  std::array<double, 4> fluxDamping{1, 1, 1, 1};
  std::array<double, 4> fluxNearFaultDamping{1, 1, 1, 1};
  // the configuration the cells of the run compute in, unless their mesh group has its own
  ConfigId config{};
  // the mesh groups whose cells compute in a configuration of their own
  std::unordered_map<int, ConfigId> groupConfigs;
  // the number of simulations of the run; the fused simulations of the cells hold them, and as
  // many copies of the cells as needed
  std::size_t simulationCount{1};

  /// The configuration the cells of a mesh group compute in.
  [[nodiscard]] ConfigId configOfGroup(int group) const;
  /// Every configuration cells of the run may compute in: `config` first, then those of the mesh
  /// groups in the order of their ids, each once.
  [[nodiscard]] std::vector<ConfigId> configs() const;
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
ModelParameters readModelParameters(ParameterReader* baseReader,
                                    ConfigId config,
                                    std::unordered_map<int, ConfigId> groupConfigs);
ITMParameters readITMParameters(ParameterReader* baseReader);
} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERS_MODELPARAMETERS_H_
