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

inline CubeGeneratorParameters readCubeGeneratorParameters(ParameterReader& baseReader) {
  auto reader = baseReader.readSubNode("cubegenerator");
  const unsigned int cubeMinX = reader.readWithDefault("cubeminx", 6);
  const unsigned int cubeMaxX = reader.readWithDefault("cubemaxx", 6);
  const unsigned int cubeMinY = reader.readWithDefault("cubeminy", 6);
  const unsigned int cubeMaxY = reader.readWithDefault("cubemaxy", 6);
  const unsigned int cubeMinZ = reader.readWithDefault("cubeminz", 6);
  const unsigned int cubeMaxZ = reader.readWithDefault("cubemaxz", 6);

  const unsigned int cubeX = reader.readWithDefault("cubex", 2);
  const unsigned int cubeY = reader.readWithDefault("cubey", 2);
  const unsigned int cubeZ = reader.readWithDefault("cubez", 2);

  // only x dimension has its number of partitions set to number of MPI processes
  const unsigned int cubePx = seissol::MPI::mpi.size();
  const unsigned int cubePy = 1;
  const unsigned int cubePz = 1;

  const double cubeS = reader.readWithDefault("cubes", 100);
  const double cubeSx = reader.readWithDefault("cubesx", cubeS);
  const double cubeSy = reader.readWithDefault("cubesy", cubeS);
  const double cubeSz = reader.readWithDefault("cubesz", cubeS);

  const double cubeTx = reader.readWithDefault("cubetx", 0.0);
  const double cubeTy = reader.readWithDefault("cubety", 0.0);
  const double cubeTz = reader.readWithDefault("cubetz", 0.0);
  return CubeGeneratorParameters{cubeMinX,
                                 cubeMaxX,
                                 cubeMinY,
                                 cubeMaxY,
                                 cubeMinZ,
                                 cubeMaxZ,
                                 cubeX,
                                 cubeY,
                                 cubeZ,
                                 cubePx,
                                 cubePy,
                                 cubePz,
                                 cubeS,
                                 cubeSx,
                                 cubeSy,
                                 cubeSz,
                                 cubeTx,
                                 cubeTy,
                                 cubeTz};
}

inline void discardCubeGeneratorParameters(ParameterReader& baseReader) {
  auto reader = baseReader.readSubNode("cubegenerator");
  reader.markUnused("cubegenerator",
                    "cubeminx",
                    "cubemaxx",
                    "cubeminy",
                    "cubemaxy",
                    "cubeminz",
                    "cubemaxz",
                    "cubex",
                    "cubey",
                    "cubez",
                    "cubes",
                    "cubesx",
                    "cubesy",
                    "cubesz",
                    "cubetx",
                    "cubety",
                    "cubetz");
}
} // namespace seissol::initializers::parameters

#endif
