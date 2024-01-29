#ifndef SEISSOL_PARAMETERS_H
#define SEISSOL_PARAMETERS_H

#include "ParameterReader.h"

#include "CubeGeneratorParameters.h"
#include "DRParameters.h"
#include "InitializationParameters.h"
#include "LtsParameters.h"
#include "MeshParameters.h"
#include "ModelParameters.h"
#include "OutputParameters.h"
#include "SourceParameters.h"

namespace seissol::initializer::parameters {

struct SeisSolParameters {
  CubeGeneratorParameters cubeGenerator;
  DRParameters drParameters;
  InitializationParameters initialization;
  MeshParameters mesh;
  ModelParameters model;
  OutputParameters output;
  SourceParameters source;
  TimeSteppingParameters timeStepping;
};

SeisSolParameters readSeisSolParameters(ParameterReader* parameterReader);
} // namespace seissol::initializer::parameters

#endif
