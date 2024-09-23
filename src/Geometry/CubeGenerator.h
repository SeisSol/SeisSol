// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 */

#ifndef SEISSOL_SRC_GEOMETRY_CUBEGENERATOR_H_
#define SEISSOL_SRC_GEOMETRY_CUBEGENERATOR_H_

#include <algorithm>
#include <cstring>
#include <map>
#include <utility>
#include <vector>

#include <omp.h>

#include "utils/logger.h"

#include "Initializer/Parameters/CubeGeneratorParameters.h"
#include "MeshReader.h"

namespace seissol::geometry {

class CubeGenerator : public seissol::geometry::MeshReader {
  int rank;
  int nProcs;

  public:
  CubeGenerator(int rank,
                int nProcs,
                const std::string& meshFile,
                const seissol::initializer::parameters::CubeGeneratorParameters& cubeParams);

  void cubeGenerator(const std::array<unsigned int, 4> numCubes,
                     const std::array<unsigned int, 4> numPartitions,
                     unsigned int boundaryMinx,
                     unsigned int boundaryMaxx,
                     unsigned int boundaryMiny,
                     unsigned int boundaryMaxy,
                     unsigned int boundaryMinz,
                     unsigned int boundaryMaxz,
                     const std::array<unsigned int, 4> numCubesPerPart,
                     const std::array<unsigned long, 4> numElemPerPart,
                     const std::array<unsigned int, 4> numVrtxPerPart,
                     const std::array<unsigned int, 3> numBndElements,
                     double scale,
                     double scaleX,
                     double scaleY,
                     double scaleZ,
                     double tx,
                     double ty,
                     double tz,
                     const std::string& meshFile);

  private:
  void findElementsPerVertex();
  void addMPINeighbor(int localID, int bndRank, int elemSize, const int* bndElemLocalIds);
};
} // namespace seissol::geometry

#endif // SEISSOL_SRC_GEOMETRY_CUBEGENERATOR_H_
