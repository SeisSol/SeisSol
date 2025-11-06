// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_GEOMETRY_CUBEGENERATOR_H_
#define SEISSOL_SRC_GEOMETRY_CUBEGENERATOR_H_

#include "Initializer/Parameters/CubeGeneratorParameters.h"
#include "MeshReader.h"
#include "utils/logger.h"

#include <algorithm>
#include <cstring>
#include <map>
#include <utility>
#include <vector>

namespace seissol::geometry {

class CubeGenerator : public seissol::geometry::MeshReader {
  int rank;
  int nProcs;

  public:
  CubeGenerator(int rank,
                int nProcs,
                const std::string& meshFile,
                const seissol::initializer::parameters::CubeGeneratorParameters& cubeParams);

  void cubeGenerator(std::array<std::size_t, 4> numCubes,
                     std::array<std::size_t, 4> numPartitions,
                     std::size_t boundaryMinx,
                     std::size_t boundaryMaxx,
                     std::size_t boundaryMiny,
                     std::size_t boundaryMaxy,
                     std::size_t boundaryMinz,
                     std::size_t boundaryMaxz,
                     std::array<std::size_t, 4> numCubesPerPart,
                     std::array<unsigned long, 4> numElemPerPart,
                     std::array<std::size_t, 4> numVrtxPerPart,
                     std::array<std::size_t, 3> numBndElements,
                     double scale,
                     double scaleX,
                     double scaleY,
                     double scaleZ,
                     double tx,
                     double ty,
                     double tz,
                     const std::string& meshFile);

  bool inlineTimestepCompute() const override;
  bool inlineClusterCompute() const override;

  private:
  void findElementsPerVertex();
  void addMPINeighbor(int localID, int bndRank, int elemSize, const int* bndElemLocalIds);
};
} // namespace seissol::geometry

#endif // SEISSOL_SRC_GEOMETRY_CUBEGENERATOR_H_
