// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "InitializationParameters.h"

#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/InputAux.h"
#include "Features/InitialConditions/InitializationParameters.h"
#include "Simulation/ConfigLayout.h"
#include "Simulation/Parameters/ModelParameters.h"
#include "Simulation/Parameters/ParameterReader.h"

#include <Eigen/Core>
#include <cstddef>
#include <string>
#include <vector>

namespace seissol::initializer::parameters {

InitializationParameters readInitializationParameters(ParameterReader* baseReader,
                                                      ConfigId config) {
  const InitializationParameters defaults{};
  auto* reader = baseReader->readSubNode("inicondition");

  const auto type = reader->readWithDefaultStringEnum<InitializationType>(
      "cictype",
      defaults.type,
      {
          {"zero", InitializationType::Zero},
          {"planarwave", InitializationType::Planarwave},
          {"superimposedplanarwave", InitializationType::SuperimposedPlanarwave},
          {"travelling", InitializationType::Travelling},
          {"acoustictravellingwithitm", InitializationType::AcousticTravellingWithITM},
          {"scholte", InitializationType::Scholte},
          {"snell", InitializationType::Snell},
          {"ocean_0", InitializationType::Ocean0},
          {"ocean_1", InitializationType::Ocean1},
          {"ocean_2", InitializationType::Ocean2},
          {"pressureinjection", InitializationType::PressureInjection},
          {"easi", InitializationType::Easi},
      });
  const auto toVector = [](const std::string& text) {
    const auto raw = seissol::initializer::convertStringToArray<double, 3>(text);
    return Eigen::Vector3d(raw.data());
  };
  const Eigen::Vector3d origin =
      reader->readWithDefaultConverted("origin", defaults.origin, toVector);
  const Eigen::Vector3d kVec = reader->readWithDefaultConverted("kvec", defaults.kVec, toVector);
  const auto numQuantities = configLayout(config).numQuantities;
  std::string defaultAmpFieldString;
  for (std::size_t i = 0; i < numQuantities; ++i) {
    defaultAmpFieldString += " 0.0";
  }
  const auto ampFieldString = reader->readWithDefault("ampfield", defaultAmpFieldString);
  const auto ampField =
      seissol::initializer::convertStringToVector<double>(ampFieldString, numQuantities);

  const auto magnitude = reader->readWithDefault("magnitude", defaults.magnitude);
  const auto width = reader->readWithDefault("width", defaults.width);
  const auto k = reader->readWithDefault("k", defaults.k);

  std::vector<std::string> filenames(readSimulationCount(baseReader, config));
  filenames[0] = reader->readPath("filename").value_or("");
  for (std::size_t i = 1; i < filenames.size(); ++i) {
    filenames[i] = reader->readPath("filename" + std::to_string(i)).value_or(filenames[0]);
  }
  const auto hasTime = reader->readWithDefault("hastime", defaults.hasTime);

  const auto avoidIC = reader->readWithDefault("avoidic", defaults.avoidIC);

  return InitializationParameters{.type = type,
                                  .origin = origin,
                                  .kVec = kVec,
                                  .ampField = ampField,
                                  .magnitude = magnitude,
                                  .width = width,
                                  .k = k,
                                  .filenames = filenames,
                                  .hasTime = hasTime,
                                  .avoidIC = avoidIC};
}
} // namespace seissol::initializer::parameters
