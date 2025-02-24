// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERS_CUBEGENERATORPARAMETERS_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERS_CUBEGENERATORPARAMETERS_H_

#include "ParameterReader.h"

namespace seissol::initializer::parameters {

struct CubeGeneratorParameters {
  unsigned int cubeMinX;
  unsigned int cubeMaxX;
  unsigned int cubeMinY;
  unsigned int cubeMaxY;
  unsigned int cubeMinZ;
  unsigned int cubeMaxZ;
  unsigned int cubeX;
  unsigned int cubeY;
  unsigned int cubeZ;
  unsigned int cubePx;
  unsigned int cubePy;
  unsigned int cubePz;
  double cubeS;
  double cubeSx;
  double cubeSy;
  double cubeSz;
  double cubeTx;
  double cubeTy;
  double cubeTz;
};

CubeGeneratorParameters readCubeGeneratorParameters(ParameterReader* baseReader);
void discardCubeGeneratorParameters(ParameterReader* baseReader);
} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERS_CUBEGENERATORPARAMETERS_H_
