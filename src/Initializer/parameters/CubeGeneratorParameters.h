#ifndef SEISSOL_CUBEGENERATOR_PARAMETERS_H
#define SEISSOL_CUBEGENERATOR_PARAMETERS_H

#include "ParameterReader.h"

namespace seissol::initializers::parameters {

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

CubeGeneratorParameters readCubeGeneratorParameters(ParameterReader& baseReader);
void discardCubeGeneratorParameters(ParameterReader& baseReader);
} // namespace seissol::initializers::parameters

#endif
