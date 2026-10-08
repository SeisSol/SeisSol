// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "InitializationParameters.h"

#include "Common/ConfigLayout.h"
#include "Common/ConfigRegistry.h"
#include "Initializer/InputAux.h"
#include "Initializer/Parameters/ParameterReader.h"

#include <Eigen/Core>
#include <cstddef>
#include <limits>
#include <string>

namespace seissol::initializer::parameters {

InitializationParameters readInitializationParameters(ParameterReader* baseReader,
                                                      ConfigId config) {
  auto* reader = baseReader->readSubNode("inicondition");

  const auto type = reader->readWithDefaultStringEnum<InitializationType>(
      "cictype",
      "zero",
      {
          {"zero", InitializationType::Zero},
          {"planarwave", InitializationType::Planarwave},
          {"superimposedplanarwave", InitializationType::SuperimposedPlanarwave},
          {"travelling", InitializationType::Travelling},
          {"scholte", InitializationType::Scholte},
          {"snell", InitializationType::Snell},
          {"ocean_0", InitializationType::Ocean0},
          {"ocean_1", InitializationType::Ocean1},
          {"ocean_2", InitializationType::Ocean2},
          {"pressureinjection", InitializationType::PressureInjection},
          {"easi", InitializationType::Script},
          {"script", InitializationType::Script},
      });
  const auto originString = reader->readWithDefault("origin", std::string("0.0 0.0 0.0"));
  const auto originRaw = seissol::initializer::convertStringToArray<double, 3>(originString);
  const Eigen::Vector3d origin(originRaw.data());
  const auto kVecString = reader->readWithDefault("kvec", std::string("0.0 0.0 0.0"));
  const auto kVecRaw = seissol::initializer::convertStringToArray<double, 3>(kVecString);
  const Eigen::Vector3d kVec(kVecRaw.data());
  const auto numQuantities = configLayout(config).numQuantities;
  std::string defaultAmpFieldString;
  for (std::size_t i = 0; i < numQuantities; ++i) {
    defaultAmpFieldString += " 0.0";
  }
  const auto ampFieldString = reader->readWithDefault("ampfield", defaultAmpFieldString);
  const auto ampField =
      seissol::initializer::convertStringToVector<double>(ampFieldString, numQuantities);

  const auto magnitude = reader->readWithDefault("magnitude", 0.0);
  const auto width = reader->readWithDefault("width", std::numeric_limits<double>::infinity());
  const auto k = reader->readWithDefault("k", 0.0);

  const auto filename = reader->readPath("filename").value_or("");
  const auto hasTime = reader->read<bool>("hastime").value_or(true);

  const auto avoidIC = reader->read<bool>("avoidic").value_or(false);

  return InitializationParameters{
      type, origin, kVec, ampField, magnitude, width, k, filename, hasTime, avoidIC};
}
} // namespace seissol::initializer::parameters
