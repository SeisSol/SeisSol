#ifndef SEISSOL_PARAMETERS_H
#define SEISSOL_PARAMETERS_H

// #include <cstdint>
// #include <string>
// #include <array>
// #include <unordered_set>
// #include <yaml-cpp/yaml.h>
//
// #include <xdmfwriter/XdmfWriter.h>
//
// #include "Checkpoint/Backend.h"
// #include "Initializer/parameters/DRParameters.h"
// #include "Geometry/CubeGenerator.h"
// #include "Geometry/MeshReader.h"
// #include "Initializer/time_stepping/LtsParameters.h"
// #include "SourceTerm/typedefs.hpp"

#include "ParameterReader.h"

#include "CubeGeneratorParameters.h"
#include "DRParameters.h"
#include "InitializationParameters.h"
#include "LtsParameters.h"
#include "MeshParameters.h"
#include "ModelParameters.h"
#include "OutputParameters.h"
#include "SourceParameters.h"

namespace seissol::initializers::parameters {

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

SeisSolParameters readSeisSolParameters(ParameterReader& parameterReader);
} // namespace seissol::initializers::parameters

#endif
