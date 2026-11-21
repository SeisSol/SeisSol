// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SIMULATION_PARAMETERS_SEISSOLPARAMETERS_H_
#define SEISSOL_SRC_SIMULATION_PARAMETERS_SEISSOLPARAMETERS_H_

#include "Core/Geometry/MeshParameters.h"
#include "Core/TimeStepping/LtsParameters.h"
#include "CubeGeneratorParameters.h"
#include "Features/DynamicRupture/DRParameters.h"
#include "Features/InitialConditions/InitializationParameters.h"
#include "Features/SourceTerm/SourceParameters.h"
#include "ModelParameters.h"
#include "OutputParameters.h"
#include "ParameterReader.h"

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

#endif // SEISSOL_SRC_SIMULATION_PARAMETERS_SEISSOLPARAMETERS_H_
