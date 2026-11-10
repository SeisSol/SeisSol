// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "SeisSolParameters.h"

#include "Common/ConfigRegistry.h"
#include "Initializer/Parameters/CubeGeneratorParameters.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Initializer/Parameters/InitializationParameters.h"
#include "Initializer/Parameters/LtsParameters.h"
#include "Initializer/Parameters/MeshParameters.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/Parameters/OutputParameters.h"
#include "Initializer/Parameters/ParameterReader.h"
#include "Initializer/Parameters/SourceParameters.h"

#include <cmath>
#include <utils/logger.h>

namespace seissol::initializer::parameters {

SeisSolParameters readSeisSolParameters(ParameterReader* parameterReader) {
  logInfo() << "Reading SeisSol parameter file...";

  // the configuration the cells of the run compute in, and the ones of the mesh groups with their
  // own; the parameters that depend on it, e.g. on its material or on its number of fused
  // simulations, are read for it
  const ConfigId config = readConfig(parameterReader);
  const auto groupConfigs = readGroupConfigs(parameterReader, config);

  const CubeGeneratorParameters cubeGeneratorParameters =
      readCubeGeneratorParameters(parameterReader);
  const DRParameters drParameters = readDRParameters(parameterReader, config);
  const InitializationParameters initializationParameters =
      readInitializationParameters(parameterReader, config);
  const MeshParameters meshParameters = readMeshParameters(parameterReader);
  const ModelParameters modelParameters =
      readModelParameters(parameterReader, config, groupConfigs);
  const OutputParameters outputParameters =
      readOutputParameters(parameterReader, modelParameters.configs());
  const SourceParameters sourceParameters = readSourceParameters(parameterReader);
  const TimeSteppingParameters timeSteppingParameters =
      readTimeSteppingParameters(parameterReader, modelParameters.configs());

  // (see physics::Ocean)
  const bool ocean = initializationParameters.type == InitializationType::Ocean0 ||
                     initializationParameters.type == InitializationType::Ocean1 ||
                     initializationParameters.type == InitializationType::Ocean2;
  if (ocean && std::abs(modelParameters.gravitationalAcceleration - 9.81e-3) > 10e-15) {
    logError() << "Ocean scenario only supports g=9.81e-3 currently!";
  }

  parameterReader->warnDeprecated({"boundaries",
                                   "rffile",
                                   "inflowbound",
                                   "inflowboundpwfile",
                                   "inflowbounduin",
                                   "source110",
                                   "source15",
                                   "source1618",
                                   "source17",
                                   "source19",
                                   "spongelayer",
                                   "sponges",
                                   "analysis",
                                   "analysisfields",
                                   "debugging"});

  logInfo() << "SeisSol parameter file read successfully.";

  return SeisSolParameters{.cubeGenerator = cubeGeneratorParameters,
                           .drParameters = drParameters,
                           .initialization = initializationParameters,
                           .mesh = meshParameters,
                           .model = modelParameters,
                           .output = outputParameters,
                           .source = sourceParameters,
                           .timeStepping = timeSteppingParameters};
}
} // namespace seissol::initializer::parameters
