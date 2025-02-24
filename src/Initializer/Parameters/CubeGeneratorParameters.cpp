// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "CubeGeneratorParameters.h"
#include <Initializer/Parameters/ParameterReader.h>

namespace seissol::initializer::parameters {

CubeGeneratorParameters readCubeGeneratorParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("cubegenerator");
  const unsigned int cubeMinX = reader->readWithDefault("cubeminx", 6);
  const unsigned int cubeMaxX = reader->readWithDefault("cubemaxx", 6);
  const unsigned int cubeMinY = reader->readWithDefault("cubeminy", 6);
  const unsigned int cubeMaxY = reader->readWithDefault("cubemaxy", 6);
  const unsigned int cubeMinZ = reader->readWithDefault("cubeminz", 6);
  const unsigned int cubeMaxZ = reader->readWithDefault("cubemaxz", 6);

  const unsigned int cubeX = reader->readWithDefault("cubex", 2);
  const unsigned int cubeY = reader->readWithDefault("cubey", 2);
  const unsigned int cubeZ = reader->readWithDefault("cubez", 2);

  // only x dimension has its number of partitions set to number of MPI processes
  const unsigned int cubePx = seissol::MPI::mpi.size();
  const unsigned int cubePy = 1;
  const unsigned int cubePz = 1;

  const double cubeS = reader->readWithDefault("cubes", 100);
  const double cubeSx = reader->readWithDefault("cubesx", cubeS);
  const double cubeSy = reader->readWithDefault("cubesy", cubeS);
  const double cubeSz = reader->readWithDefault("cubesz", cubeS);

  const double cubeTx = reader->readWithDefault("cubetx", 0.0);
  const double cubeTy = reader->readWithDefault("cubety", 0.0);
  const double cubeTz = reader->readWithDefault("cubetz", 0.0);
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

void discardCubeGeneratorParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("cubegenerator");
  reader->markUnused({"cubegenerator",
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
                      "cubetz"});
}
} // namespace seissol::initializer::parameters
