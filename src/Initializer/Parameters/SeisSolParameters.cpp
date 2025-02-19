// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "SeisSolParameters.h"
#include <Initializer/Parameters/CubeGeneratorParameters.h>
#include <Initializer/Parameters/DRParameters.h>
#include <Initializer/Parameters/InitializationParameters.h>
#include <Initializer/Parameters/LtsParameters.h>
#include <Initializer/Parameters/MeshParameters.h>
#include <Initializer/Parameters/ModelParameters.h>
#include <Initializer/Parameters/OutputParameters.h>
#include <Initializer/Parameters/ParameterReader.h>
#include <Initializer/Parameters/SourceParameters.h>
#include <utils/logger.h>

namespace seissol::initializer::parameters {

SeisSolParameters readSeisSolParameters(ParameterReader* parameterReader) {
  logInfo() << "Reading SeisSol parameter file...";

  const CubeGeneratorParameters cubeGeneratorParameters =
      readCubeGeneratorParameters(parameterReader);
  const DRParameters drParameters = readDRParameters(parameterReader);
  const InitializationParameters initializationParameters =
      readInitializationParameters(parameterReader);
  const MeshParameters meshParameters = readMeshParameters(parameterReader);
  const ModelParameters modelParameters = readModelParameters(parameterReader);
  const OutputParameters outputParameters = readOutputParameters(parameterReader);
  const SourceParameters sourceParameters = readSourceParameters(parameterReader);
  const TimeSteppingParameters timeSteppingParameters = readTimeSteppingParameters(parameterReader);

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

  return SeisSolParameters{cubeGeneratorParameters,
                           drParameters,
                           initializationParameters,
                           meshParameters,
                           modelParameters,
                           outputParameters,
                           sourceParameters,
                           timeSteppingParameters};
}
} // namespace seissol::initializer::parameters
