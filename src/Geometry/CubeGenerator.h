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

#include "utils/logger.h"
#include "utils/args.h"
#include "MeshReader.h"

#include <omp.h>

#include <algorithm>
#include <cstring>
#include <map>
#include <utility>
#include <vector>

#include <iostream>

namespace seissol::geometry {

struct CubeGeneratorParameters {
  unsigned int cubeMinX;
  unsigned int cubeMaxX;
  unsigned int cubeMinY;
  unsigned int cubeMaxY;
  unsigned int cubeMinZ;
  unsigned int cubeMaxZ;
  unsigned cubeX;
  unsigned cubeY;
  unsigned cubeZ;
  unsigned cubePx;
  unsigned cubePy;
  unsigned cubePz;
  double cubeS;
  double cubeSx;
  double cubeSy;
  double cubeSz;
  double cubeTx;
  double cubeTy;
  double cubeTz;
};

class CubeGenerator : public seissol::geometry::MeshReader {
  public:
  CubeGenerator(int rank,
                int nProcs,
                const std::string& meshFile,
                const seissol::geometry::CubeGeneratorParameters& cubeParams);

  /*
    inline void loadBar(int x, int n, int r = 100, int w = 50);
    const char* dim2str(unsigned int dim);
    template <typename A, typename B>
    std::pair<B, A> flip_pair(const std::pair<A, B>& p);
    void checkNcError(int error);
  */
  void cubeGenerator(unsigned int numCubes[4],
                     unsigned int numPartitions[4],
                     unsigned int boundaryMinx,
                     unsigned int boundaryMaxx,
                     unsigned int boundaryMiny,
                     unsigned int boundaryMaxy,
                     unsigned int boundaryMinz,
                     unsigned int boundaryMaxz,
                     unsigned int numCubesPerPart[4],
                     unsigned long numElemPerPart[4],
                     unsigned int numVrtxPerPart[4],
                     unsigned int numBndElements[3],
                     double scale,
                     double scaleX,
                     double scaleY,
                     double scaleZ,
                     double tx,
                     double ty,
                     double tz,
                     const std::string& output);
};
} // namespace seissol::geometry
#endif // CUBEGENERATOR_H
