/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
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
 */

#ifndef PUMLREADER_H
#define PUMLREADER_H

#include "Initializer/Parameters/MeshParameters.h"
#include "MeshReader.h"
#include "PUML/PUML.h"
#include "Parallel/MPI.h"

namespace seissol::initializer::time_stepping {
class LtsWeights;
} // namespace seissol::initializer::time_stepping

namespace seissol::geometry {
inline int decodeBoundary(const void* data,
                          size_t cell,
                          int face,
                          seissol::initializer::parameters::BoundaryFormat format) {
  if (format == seissol::initializer::parameters::BoundaryFormat::I32) {
    const auto* dataCasted = reinterpret_cast<const uint32_t*>(data);
    return (dataCasted[cell] >> (8 * face)) & 0xff;
  } else if (format == seissol::initializer::parameters::BoundaryFormat::I64) {
    const auto* dataCasted = reinterpret_cast<const uint64_t*>(data);
    return (dataCasted[cell] >> (16 * face)) & 0xffff;
  } else if (format == seissol::initializer::parameters::BoundaryFormat::I32x4) {
    const int* dataCasted = reinterpret_cast<const int*>(data);
    return dataCasted[cell * 4 + face];
  } else {
    logError() << "Unknown boundary format:" << static_cast<int>(format);
    return 0;
  }
}

class PUMLReader : public seissol::geometry::MeshReader {
  public:
  PUMLReader(const char* meshFile,
             const char* partitioningLib,
             double maximumAllowedTimeStep,
             seissol::initializer::parameters::BoundaryFormat boundaryFormat =
                 seissol::initializer::parameters::BoundaryFormat::I32,
             initializer::time_stepping::LtsWeights* ltsWeights = nullptr,
             double tpwgt = 1.0);

  private:
  seissol::initializer::parameters::BoundaryFormat boundaryFormat;

  /**
   * Read the mesh
   */
  void read(PUML::TETPUML& puml, const char* meshFile);

  /**
   * Create the partitioning
   */
  static void partition(PUML::TETPUML& puml,
                        initializer::time_stepping::LtsWeights* ltsWeights,
                        double tpwgt,
                        const char* meshFile,
                        const char* partitioningLib);
  /**
   * Generate the PUML data structure
   */
  static void generatePUML(PUML::TETPUML& puml);

  /**
   * Get the mesh
   */
  void getMesh(const PUML::TETPUML& puml);

  void addMPINeighor(const PUML::TETPUML& puml, int rank, const std::vector<unsigned int>& faces);
};

} // namespace seissol::geometry

#endif // PUMLREADER_H
