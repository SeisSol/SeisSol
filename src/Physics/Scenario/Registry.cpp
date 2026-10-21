// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Registry.h"

#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"
#include "Initializer/Parameters/InitializationParameters.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/Typedefs.h"
#include "Model/CommonDatastructures.h"
#include "Model/MaterialType.h"
#include "Physics/InitialField.h"
#include "Physics/Scenario/PlaneWave.h"
#include "Physics/Scenario/Scenarios.h"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <math.h>
#include <memory>
#include <string_view>
#include <vector>

namespace seissol::physics::scenario {

namespace {

using initializer::parameters::InitializationType;

using FieldList = std::vector<std::unique_ptr<InitialField>>;

// The material models a scenario is formulated for.
enum class MaterialRequirement { None, NotAnelastic, Poroelastic, Fluid, SolidAndFluid };

bool fulfills(const ConfigValue& config, MaterialRequirement requirement) {
  switch (requirement) {
  case MaterialRequirement::NotAnelastic:
    return config.relaxationMechanisms == 0;
  case MaterialRequirement::Poroelastic:
    return config.materialType == model::MaterialType::Poroelastic;
  case MaterialRequirement::Fluid:
    // the cells with mu = 0 of the elastic material, or the acoustic one; the closed-form
    // solutions are lossless, so the attenuating materials are left out
    return config.materialType == model::MaterialType::Elastic ||
           config.materialType == model::MaterialType::Acoustic;
  case MaterialRequirement::SolidAndFluid:
    return config.materialType == model::MaterialType::Elastic;
  case MaterialRequirement::None:
    break;
  }
  return true;
}

TravellingWaveParameters travellingWaveParameters(const Input& input) {
  const auto& initialization = input.parameters.initialization;

  TravellingWaveParameters parameters{};
  parameters.origin = initialization.origin;
  parameters.kVec = initialization.kVec;
  constexpr double Eps = 1e-15;
  for (std::size_t i = 0; i < initialization.ampField.size(); ++i) {
    if (std::abs(initialization.ampField[i]) > Eps) {
      parameters.varField.push_back(i);
      parameters.ampField.emplace_back(initialization.ampField[i]);
    }
  }
  return parameters;
}

struct Entry {
  InitializationType type;
  std::string_view name;
  // What the material has to be for the scenario to be defined.
  MaterialRequirement requirement;
  // Full clause, used as-is in the diagnostic when the scenario is unavailable.
  std::string_view materialRequirement;
  // Whether the scenario can be evaluated at arbitrary times, which an analytical boundary
  // condition needs.
  bool timeDependent;
  FieldList (*build)(const Input&);
};

constexpr std::array<Entry, 12> Registry = {{
    {InitializationType::Zero,
     "zero",
     MaterialRequirement::None,
     "",
     true,
     [](const Input& input) { return zeroField(input.config); }},
    {InitializationType::Planarwave,
     "planar wave",
     MaterialRequirement::None,
     "",
     true,
     // along kVec where the parameter file sets it, else along (pi, pi, pi); fused simulations
     // are offset against each other in phase, so that they do not all carry the same wave
     [](const Input& input) {
       return planarwave(input.materialData, input.config, input.parameters.initialization.kVec);
     }},
    {InitializationType::SuperimposedPlanarwave,
     "super-imposed planar wave",
     MaterialRequirement::None,
     "",
     true,
     [](const Input& input) { return superimposedPlanarwave(input.materialData, input.config); }},
    {InitializationType::Travelling,
     "travelling wave",
     MaterialRequirement::NotAnelastic,
     "it is formulated for a material without anelastic mechanisms",
     true,
     [](const Input& input) {
       return travellingWave(input.materialData, input.config, travellingWaveParameters(input));
     }},
    {InitializationType::AcousticTravellingWithITM,
     "acoustic travelling wave with ITM",
     MaterialRequirement::Fluid,
     "it is formulated for a fluid, i.e. the elastic material with mu = 0 or the acoustic one",
     true,
     [](const Input& input) {
       return acousticTravellingWaveITM(
           input.materialData,
           acousticTravellingWaveParameters(input.parameters.initialization,
                                            input.parameters.model.itmParameters),
           input.config);
     }},
    {InitializationType::Scholte,
     "Scholte wave (elastic-acoustic)",
     MaterialRequirement::SolidAndFluid,
     "it couples a solid to a fluid, which only the elastic material holds both of",
     true,
     [](const Input& input) { return scholteWave(input.config); }},
    {InitializationType::Snell,
     "Snell's law (elastic-acoustic)",
     MaterialRequirement::SolidAndFluid,
     "it couples a solid to a fluid, which only the elastic material holds both of",
     true,
     [](const Input& input) { return snellsLaw(input.config); }},
    {InitializationType::Ocean0,
     "ocean, an uncoupled ocean test case for acoustic equations (mode 0)",
     MaterialRequirement::Fluid,
     "it is formulated for a fluid, i.e. the elastic material with mu = 0 or the acoustic one",
     true,
     [](const Input& input) { return ocean(0, input.gravitation.acceleration, input.config); }},
    {InitializationType::Ocean1,
     "ocean, an uncoupled ocean test case for acoustic equations (mode 1)",
     MaterialRequirement::Fluid,
     "it is formulated for a fluid, i.e. the elastic material with mu = 0 or the acoustic one",
     true,
     [](const Input& input) { return ocean(1, input.gravitation.acceleration, input.config); }},
    {InitializationType::Ocean2,
     "ocean, an uncoupled ocean test case for acoustic equations (mode 2)",
     MaterialRequirement::Fluid,
     "it is formulated for a fluid, i.e. the elastic material with mu = 0 or the acoustic one",
     true,
     [](const Input& input) { return ocean(2, input.gravitation.acceleration, input.config); }},
    {InitializationType::PressureInjection,
     "pressure injection",
     MaterialRequirement::Poroelastic,
     "it is formulated for a poroelastic material",
     true,
     [](const Input& input) {
       return pressureInjection(input.parameters.initialization, input.config);
     }},
    // Read from a script, which may read the time t: physics::ScriptField evaluates it at any
    // point and time, as an analytical boundary condition asks.
    {InitializationType::Script, "script", MaterialRequirement::None, {}, true, nullptr},
}};

const Entry* find(InitializationType type) {
  const auto* entry = std::find_if(
      Registry.begin(), Registry.end(), [type](const Entry& e) { return e.type == type; });
  return entry == Registry.end() ? nullptr : entry;
}

} // namespace

std::string_view name(InitializationType type) {
  const auto* entry = find(type);
  return entry == nullptr ? "unknown" : entry->name;
}

Availability availability(InitializationType type, ConfigId config) {
  const auto* entry = find(type);
  if (entry == nullptr) {
    return {false, "it is not a known scenario"};
  }
  if (!fulfills(configValue(config), entry->requirement)) {
    return {false, entry->materialRequirement};
  }
  return {};
}

Availability analyticalBoundaryAvailability(InitializationType type, ConfigId config) {
  const auto general = availability(type, config);
  if (!general.available) {
    return general;
  }
  const auto* entry = find(type);
  if (!entry->timeDependent) {
    return {false, "it supplies a state at t = 0 only"};
  }
  return {};
}

AcousticTravellingWaveParametersITM acousticTravellingWaveParameters(
    const initializer::parameters::InitializationParameters& initialization,
    const initializer::parameters::ITMParameters& itm) {
  AcousticTravellingWaveParametersITM parameters{};
  parameters.k = initialization.k;
  parameters.itmStartingTime = itm.itmStartingTime;
  parameters.itmDuration = itm.itmDuration;
  // Without an enabled mirror, the ITM values are ignored. A factor of one makes the wave pass the
  // mirror times unchanged, whatever they are.
  parameters.itmVelocityScalingFactor = itm.itmEnabled ? itm.itmVelocityScalingFactor : 1.0;
  return parameters;
}

FieldList build(InitializationType type, const Input& input) {
  const auto* entry = find(type);
  if (entry == nullptr || entry->build == nullptr) {
    return {};
  }
  return entry->build(input);
}

} // namespace seissol::physics::scenario
