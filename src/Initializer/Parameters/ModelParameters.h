// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERS_MODELPARAMETERS_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERS_MODELPARAMETERS_H_

#include "Common/ConfigRegistry.h"
#include "Common/SegmentMap.h"
#include "ParameterReader.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
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

//! @brief The flux with its parts (see ModelParameters::fluxGPart and the following), as the log
//! shows it.
std::string describeFlux(NumericalFlux flux,
                         double gpart,
                         double rpart,
                         const std::array<double, 6>& damping);

/**
 * A flux with its parts: the share of the Godunov flux (the rest is the centered flux), that of the
 * damping of the Rusanov flux, which only "mixed" and its presets have, and the shares of the
 * upwind damping in the Godunov flux (see model::DampingShares).
 */
struct FluxSetting {
  NumericalFlux flux{NumericalFlux::Godunov};
  double gpart{0};
  double rpart{0};
  std::array<double, 6> damping{1, 1, 1, 1, 1, 1};

  bool operator==(const FluxSetting& other) const = default;

  [[nodiscard]] std::string describe() const { return describeFlux(flux, gpart, rpart, damping); }
};

//! @brief The flux that numflux<suffix> names in @p reader , with its parts fluxgpart<suffix>,
//! fluxrpart<suffix> and fluxdamping<suffix>; the parameter file and the entries of the face map
//! give them alike (suffix "" or "nearfault").
FluxSetting readFluxSetting(ParameterReader* reader, const std::string& suffix);

/**
 * A rule of the face map for the flux of regular faces: it selects the faces with one of its tags
 * (any, without tags) that have a cell of one of its mesh groups on either side (any, without
 * groups), and sets the flux of those near the fault, of the others, or of both.
 */
struct FaceFluxRule {
  std::optional<SegmentMap<std::uint32_t, bool>> tags;
  std::optional<SegmentMap<int, bool>> groups;
  std::optional<FluxSetting> flux;
  std::optional<FluxSetting> fluxNearFault;
  // the entry of the face map, for the log
  std::string name;

  [[nodiscard]] bool selects(std::uint32_t tag, int group, int neighborGroup) const;
};

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
  // the shares of the upwind damping in the Godunov part of the flux (see model::DampingShares):
  // of the velocity jumps (normal, tangential), of the traction jumps (normal, tangential), and of
  // the jumps of the velocity and the pressure of a fluid (normal); all 1 is the Godunov flux
  std::array<double, 6> fluxDamping{1, 1, 1, 1, 1, 1};
  std::array<double, 6> fluxNearFaultDamping{1, 1, 1, 1, 1, 1};
  // the rules of the face map for the flux of regular faces, in its order (see FaceFluxRule)
  std::vector<FaceFluxRule> fluxRules;
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
  /// The flux the parameter file gives the regular faces near the fault, or the others.
  [[nodiscard]] FluxSetting fluxSetting(bool nearFault) const;
  /// The rule of the face map that decides the flux of a regular face with the tag `tag` between
  /// cells of the mesh groups `group` and `neighborGroup`: the first that selects the face and sets
  /// a flux for it (near the fault, or not); none if the parameter file does.
  [[nodiscard]] std::optional<std::size_t>
      fluxRuleOf(std::uint32_t tag, int group, int neighborGroup, bool nearFault) const;
  /// The flux of such a face, by its rule or by the parameter file.
  [[nodiscard]] FluxSetting
      fluxOf(std::uint32_t tag, int group, int neighborGroup, bool nearFault) const;
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
