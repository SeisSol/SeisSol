/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 **/
#ifndef CUBEGENERATOR_H
#define CUBEGENERATOR_H

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

  void cubeGenerator(std::array<unsigned int, 4> numCubes,
                     std::array<unsigned int, 4> numPartitions,
                     unsigned int boundaryMinx,
                     unsigned int boundaryMaxx,
                     unsigned int boundaryMiny,
                     unsigned int boundaryMaxy,
                     unsigned int boundaryMinz,
                     unsigned int boundaryMaxz,
                     std::array<unsigned int, 4> numCubesPerPart,
                     std::array<unsigned long, 4> numElemPerPart,
                     std::array<unsigned int, 4> numVrtxPerPart,
                     std::array<unsigned int, 3> numBndElements,
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
#endif // CUBEGENERATOR_H
