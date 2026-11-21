// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "SeisSolParameters.h"

#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/RunConfigs.h"
#include "Core/Geometry/MeshParameters.h"
#include "Core/TimeStepping/LtsParameters.h"
#include "Features/DynamicRupture/DRParameters.h"
#include "Features/InitialConditions/InitializationParameters.h"
#include "Features/Plasticity/PlasticityParameters.h"
#include "Features/SourceTerm/SourceParameters.h"
#include "Features/TimeMirror/ITMParameters.h"
#include "Simulation/Parameters/CubeGeneratorParameters.h"
#include "Simulation/Parameters/DRParameters.h"
#include "Simulation/Parameters/InitializationParameters.h"
#include "Simulation/Parameters/LtsParameters.h"
#include "Simulation/Parameters/MeshParameters.h"
#include "Simulation/Parameters/ModelParameters.h"
#include "Simulation/Parameters/OutputParameters.h"
#include "Simulation/Parameters/ParameterReader.h"
#include "Simulation/Parameters/SourceParameters.h"
#include "WavePropagation/WavePropagationParameters.h"

#include <cmath>
#include <utils/logger.h>

namespace seissol::initializer::parameters {

SeisSolParameters readSeisSolParameters(ParameterReader* parameterReader) {
  logInfo() << "Reading SeisSol parameter file...";

  // the configuration the cells of the run compute in, and the ones of the mesh groups with their
  // own; the parameters that depend on it, e.g. on its material or on its number of fused
  // simulations, are read for it
  const RunConfigs runConfigs = readRunConfigs(parameterReader);
  const ConfigId config = runConfigs.config;

  const CubeGeneratorParameters cubeGeneratorParameters =
      readCubeGeneratorParameters(parameterReader);
  const DRParameters drParameters = readDRParameters(parameterReader, config);
  const InitializationParameters initializationParameters =
      readInitializationParameters(parameterReader, config);
  const MeshParameters meshParameters = readMeshParameters(parameterReader);
  const ModelParameters modelParameters = readModelParameters(parameterReader, runConfigs);
  const WavePropagationParameters wavePropagationParameters =
      readWavePropagationParameters(parameterReader, runConfigs);
  const PlasticityParameters plasticityParameters = readPlasticityParameters(parameterReader);
  const ITMParameters timeMirrorParameters = readITMParameters(parameterReader);
  const OutputParameters outputParameters = readOutputParameters(
      parameterReader, runConfigs.configs(), drParameters.isDynamicRuptureEnabled);
  const SourceParameters sourceParameters = readSourceParameters(parameterReader);
  const TimeSteppingParameters timeSteppingParameters =
      readTimeSteppingParameters(parameterReader, runConfigs.configs());

  // (see physics::Ocean)
  const bool ocean = initializationParameters.type == InitializationType::Ocean0 ||
                     initializationParameters.type == InitializationType::Ocean1 ||
                     initializationParameters.type == InitializationType::Ocean2;
  if (ocean && std::abs(wavePropagationParameters.gravitationalAcceleration - 9.81e-3) > 10e-15) {
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
                           .plasticity = plasticityParameters,
                           .runConfigs = runConfigs,
                           .source = sourceParameters,
                           .timeMirror = timeMirrorParameters,
                           .timeStepping = timeSteppingParameters,
                           .wavePropagation = wavePropagationParameters};
}
} // namespace seissol::initializer::parameters
