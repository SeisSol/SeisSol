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
#include "Physics/Scenario/Scenarios.h"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
// M_PI is POSIX and only declared by the C header
#include <math.h> // NOLINT(modernize-deprecated-headers)
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

// Fused simulations are offset against each other in phase, so that they do not all carry the
// same wave.
double simulationPhase(std::size_t sim, std::size_t numSimulations) {
  return (2.0 * M_PI * sim) / numSimulations;
}

template <typename FieldT, typename... Args>
FieldList perSimulation(const Input& input, const Args&... args) {
  FieldList fields;
  const auto numSimulations = configValue(input.config).numSimulations;
  for (std::size_t sim = 0; sim < numSimulations; ++sim) {
    fields.emplace_back(std::make_unique<FieldT>(args..., simulationPhase(sim, numSimulations)));
  }
  return fields;
}

// The planar wave travels along kVec where the parameter file sets it, and keeps its default wave
// vector otherwise.
FieldList planarWaves(const Input& input) {
  const Eigen::Vector3d& kVec = input.parameters.initialization.kVec;
  if (kVec.isZero(0.0)) {
    return perSimulation<Planarwave>(input, input.materialData, input.config);
  }
  FieldList fields;
  const auto numSimulations = configValue(input.config).numSimulations;
  for (std::size_t sim = 0; sim < numSimulations; ++sim) {
    fields.emplace_back(std::make_unique<Planarwave>(
        input.materialData, input.config, simulationPhase(sim, numSimulations), kVec));
  }
  return fields;
}

template <typename FieldT, typename... Args>
FieldList single(Args&&... args) {
  FieldList fields;
  fields.emplace_back(std::make_unique<FieldT>(std::forward<Args>(args)...));
  return fields;
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
    {.type = InitializationType::Zero,
     .name = "zero",
     .requirement = MaterialRequirement::None,
     .materialRequirement = "",
     .timeDependent = true,
     .build = [](const Input&) { return single<ZeroField>(); }},
    {.type = InitializationType::Planarwave,
     .name = "planar wave",
     .requirement = MaterialRequirement::None,
     .materialRequirement = "",
     .timeDependent = true,
     .build = [](const Input& input) { return planarWaves(input); }},
    {.type = InitializationType::SuperimposedPlanarwave,
     .name = "super-imposed planar wave",
     .requirement = MaterialRequirement::None,
     .materialRequirement = "",
     .timeDependent = true,
     .build =
         [](const Input& input) {
           return perSimulation<SuperimposedPlanarwave>(input, input.materialData, input.config);
         }},
    {.type = InitializationType::Travelling,
     .name = "travelling wave",
     .requirement = MaterialRequirement::NotAnelastic,
     .materialRequirement = "it is formulated for a material without anelastic mechanisms",
     .timeDependent = true,
     .build =
         [](const Input& input) {
           return single<TravellingWave>(
               input.materialData, input.config, travellingWaveParameters(input));
         }},
    {.type = InitializationType::AcousticTravellingWithITM,
     .name = "acoustic travelling wave with ITM",
     .requirement = MaterialRequirement::Fluid,
     .materialRequirement =
         "it is formulated for a fluid, i.e. the elastic material with mu = 0 or the acoustic one",
     .timeDependent = true,
     .build =
         [](const Input& input) {
           return single<AcousticTravellingWaveITM>(
               input.materialData,
               acousticTravellingWaveParameters(input.parameters.initialization,
                                                input.parameters.model.itmParameters),
               input.config);
         }},
    {.type = InitializationType::Scholte,
     .name = "Scholte wave (elastic-acoustic)",
     .requirement = MaterialRequirement::SolidAndFluid,
     .materialRequirement =
         "it couples a solid to a fluid, which only the elastic material holds both of",
     .timeDependent = true,
     .build = [](const Input&) { return single<ScholteWave>(); }},
    {.type = InitializationType::Snell,
     .name = "Snell's law (elastic-acoustic)",
     .requirement = MaterialRequirement::SolidAndFluid,
     .materialRequirement =
         "it couples a solid to a fluid, which only the elastic material holds both of",
     .timeDependent = true,
     .build = [](const Input&) { return single<SnellsLaw>(); }},
    {.type = InitializationType::Ocean0,
     .name = "ocean, an uncoupled ocean test case for acoustic equations (mode 0)",
     .requirement = MaterialRequirement::Fluid,
     .materialRequirement =
         "it is formulated for a fluid, i.e. the elastic material with mu = 0 or the acoustic one",
     .timeDependent = true,
     .build =
         [](const Input& input) {
           return single<Ocean>(0, input.gravitation.acceleration, input.config);
         }},
    {.type = InitializationType::Ocean1,
     .name = "ocean, an uncoupled ocean test case for acoustic equations (mode 1)",
     .requirement = MaterialRequirement::Fluid,
     .materialRequirement =
         "it is formulated for a fluid, i.e. the elastic material with mu = 0 or the acoustic one",
     .timeDependent = true,
     .build =
         [](const Input& input) {
           return single<Ocean>(1, input.gravitation.acceleration, input.config);
         }},
    {.type = InitializationType::Ocean2,
     .name = "ocean, an uncoupled ocean test case for acoustic equations (mode 2)",
     .requirement = MaterialRequirement::Fluid,
     .materialRequirement =
         "it is formulated for a fluid, i.e. the elastic material with mu = 0 or the acoustic one",
     .timeDependent = true,
     .build =
         [](const Input& input) {
           return single<Ocean>(2, input.gravitation.acceleration, input.config);
         }},
    {.type = InitializationType::PressureInjection,
     .name = "pressure injection",
     .requirement = MaterialRequirement::Poroelastic,
     .materialRequirement = "it is formulated for a poroelastic material",
     .timeDependent = true,
     .build =
         [](const Input& input) {
           return single<PressureInjection>(input.parameters.initialization);
         }},
    // Read from a file, and hence known at t = 0 only.
    {.type = InitializationType::Easi,
     .name = "easi file",
     .requirement = MaterialRequirement::None,
     .materialRequirement = {},
     .timeDependent = false,
     .build = nullptr},
}};

const Entry* find(InitializationType type) {
  const auto* entry =
      std::ranges::find_if(Registry, [type](const Entry& e) { return e.type == type; });
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
    return {.available = false, .reason = "it is not a known scenario"};
  }
  if (!fulfills(configValue(config), entry->requirement)) {
    return {.available = false, .reason = entry->materialRequirement};
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
    return {.available = false, .reason = "it supplies a state at t = 0 only"};
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
