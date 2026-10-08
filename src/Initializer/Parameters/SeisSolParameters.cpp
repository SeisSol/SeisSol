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
#include "Initializer/Parameters/DatafieldParameters.h"
#include "Initializer/Parameters/InitializationParameters.h"
#include "Initializer/Parameters/LtsParameters.h"
#include "Initializer/Parameters/MeshParameters.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/Parameters/OutputParameters.h"
#include "Initializer/Parameters/ParameterReader.h"
#include "Initializer/Parameters/SourceParameters.h"

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
  const DatafieldParameters datafieldParameters = readDatafieldParameters(parameterReader);
  const DRParameters drParameters = readDRParameters(parameterReader, config);
  const InitializationParameters initializationParameters =
      readInitializationParameters(parameterReader, config);
  const MeshParameters meshParameters = readMeshParameters(parameterReader);
  const ModelParameters modelParameters =
      readModelParameters(parameterReader, config, groupConfigs);
  const OutputParameters outputParameters = readOutputParameters(parameterReader, config);
  const SourceParameters sourceParameters = readSourceParameters(parameterReader);
  const TimeSteppingParameters timeSteppingParameters =
      readTimeSteppingParameters(parameterReader, modelParameters.configs());

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
                           datafieldParameters,
                           drParameters,
                           initializationParameters,
                           meshParameters,
                           modelParameters,
                           outputParameters,
                           sourceParameters,
                           timeSteppingParameters};
}
} // namespace seissol::initializer::parameters
