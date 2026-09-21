// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Registry.h"

#include "Equations/Datastructures.h"
#include "Initializer/Parameters/InitializationParameters.h"
#include "Initializer/Typedefs.h"
#include "Model/CommonDatastructures.h"
#include "Physics/InitialField.h"
#include "Physics/Scenario/Scenarios.h"
#include "Solver/MultipleSimulations.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <string_view>
#include <vector>

namespace seissol::physics::scenario {

namespace {

using initializer::parameters::InitializationType;

using FieldList = std::vector<std::unique_ptr<InitialField>>;

// The material models a scenario is formulated for. Kept as a compile-time condition rather
// than a runtime check, since the material is fixed when the binary is built.
constexpr bool Anelastic = model::MaterialT::Mechanisms > 0;
constexpr bool Poroelastic = model::MaterialT::Type == model::MaterialType::Poroelastic;

TravellingWaveParameters travellingWaveParameters(const Input& input) {
  const auto& initialization = input.parameters.initialization;

  TravellingWaveParameters parameters{};
  parameters.origin = initialization.origin;
  parameters.kVec = initialization.kVec;
  constexpr double Eps = 1e-15;
  for (std::size_t i = 0; i < model::MaterialT::NumQuantities; ++i) {
    if (std::abs(initialization.ampField[i]) > Eps) {
      parameters.varField.push_back(i);
      parameters.ampField.emplace_back(initialization.ampField[i]);
    }
  }
  return parameters;
}

AcousticTravellingWaveParametersITM acousticTravellingWaveParameters(const Input& input) {
  const auto& initialization = input.parameters.initialization;
  const auto& itm = input.parameters.model.itmParameters;

  AcousticTravellingWaveParametersITM parameters{};
  parameters.k = initialization.k;
  parameters.itmStartingTime = itm.itmStartingTime;
  parameters.itmDuration = itm.itmDuration;
  parameters.itmVelocityScalingFactor = itm.itmVelocityScalingFactor;
  return parameters;
}

// Fused simulations are offset against each other in phase, so that they do not all carry the
// same wave.
template <typename FieldT, typename... Args>
FieldList perSimulation(Args&&... args) {
  FieldList fields;
  for (std::size_t sim = 0; sim < multisim::NumSimulations; ++sim) {
    const double phase = (2.0 * M_PI * sim) / multisim::NumSimulations;
    fields.emplace_back(std::make_unique<FieldT>(args..., phase));
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
  // Whether the scenario is defined for the material this binary was built for.
  bool definedForMaterial;
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
     true,
     "",
     true,
     [](const Input&) { return single<ZeroField>(); }},
    {InitializationType::Planarwave,
     "planar wave",
     true,
     "",
     true,
     [](const Input& input) { return perSimulation<Planarwave>(input.materialData); }},
    {InitializationType::SuperimposedPlanarwave,
     "super-imposed planar wave",
     true,
     "",
     true,
     [](const Input& input) { return perSimulation<SuperimposedPlanarwave>(input.materialData); }},
    {InitializationType::Travelling,
     "travelling wave",
     !Anelastic,
     "it is formulated for a material without anelastic mechanisms",
     true,
     [](const Input& input) {
       return single<TravellingWave>(input.materialData, travellingWaveParameters(input));
     }},
    {InitializationType::AcousticTravellingWithITM,
     "acoustic travelling wave with ITM",
     !Anelastic,
     "it is formulated for a material without anelastic mechanisms",
     true,
     [](const Input& input) {
       return single<AcousticTravellingWaveITM>(input.materialData,
                                                acousticTravellingWaveParameters(input));
     }},
    {InitializationType::Scholte,
     "Scholte wave (elastic-acoustic)",
     !Anelastic,
     "it is formulated for a material without anelastic mechanisms",
     true,
     [](const Input&) { return single<ScholteWave>(); }},
    {InitializationType::Snell,
     "Snell's law (elastic-acoustic)",
     !Anelastic,
     "it is formulated for a material without anelastic mechanisms",
     true,
     [](const Input&) { return single<SnellsLaw>(); }},
    {InitializationType::Ocean0,
     "ocean, an uncoupled ocean test case for acoustic equations (mode 0)",
     !Anelastic,
     "it is formulated for a material without anelastic mechanisms",
     true,
     [](const Input& input) { return single<Ocean>(0, input.gravitation.acceleration); }},
    {InitializationType::Ocean1,
     "ocean, an uncoupled ocean test case for acoustic equations (mode 1)",
     !Anelastic,
     "it is formulated for a material without anelastic mechanisms",
     true,
     [](const Input& input) { return single<Ocean>(1, input.gravitation.acceleration); }},
    {InitializationType::Ocean2,
     "ocean, an uncoupled ocean test case for acoustic equations (mode 2)",
     !Anelastic,
     "it is formulated for a material without anelastic mechanisms",
     true,
     [](const Input& input) { return single<Ocean>(2, input.gravitation.acceleration); }},
    {InitializationType::PressureInjection,
     "pressure injection",
     Poroelastic,
     "it is formulated for a poroelastic material",
     true,
     [](const Input& input) { return single<PressureInjection>(input.parameters.initialization); }},
    // Read from a file, and hence known at t = 0 only.
    {InitializationType::Easi, "easi file", true, {}, false, nullptr},
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

Availability availability(InitializationType type) {
  const auto* entry = find(type);
  if (entry == nullptr) {
    return {false, "it is not a known scenario"};
  }
  if (!entry->definedForMaterial) {
    return {false, entry->materialRequirement};
  }
  return {};
}

Availability analyticalBoundaryAvailability(InitializationType type) {
  const auto general = availability(type);
  if (!general.available) {
    return general;
  }
  const auto* entry = find(type);
  if (!entry->timeDependent) {
    return {false, "it supplies a state at t = 0 only"};
  }
  return {};
}

FieldList build(InitializationType type, const Input& input) {
  const auto* entry = find(type);
  if (entry == nullptr || entry->build == nullptr) {
    return {};
  }
  return entry->build(input);
}

} // namespace seissol::physics::scenario
