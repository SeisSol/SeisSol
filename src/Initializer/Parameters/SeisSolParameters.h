// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERS_SEISSOLPARAMETERS_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERS_SEISSOLPARAMETERS_H_

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

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERS_SEISSOLPARAMETERS_H_
