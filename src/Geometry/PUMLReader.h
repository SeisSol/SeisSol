// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 */

#ifndef SEISSOL_SRC_GEOMETRY_PUMLREADER_H_
#define SEISSOL_SRC_GEOMETRY_PUMLREADER_H_

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
    const uint32_t* dataCasted = reinterpret_cast<const uint32_t*>(data);
    return (dataCasted[cell] >> (8 * face)) & 0xff;
  } else if (format == seissol::initializer::parameters::BoundaryFormat::I64) {
    const uint64_t* dataCasted = reinterpret_cast<const uint64_t*>(data);
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
             const char* checkPointFile,
             seissol::initializer::parameters::BoundaryFormat boundaryFormat =
                 seissol::initializer::parameters::BoundaryFormat::I32,
             initializer::time_stepping::LtsWeights* ltsWeights = nullptr,
             double tpwgt = 1.0,
             bool readPartitionFromFile = false);

  private:
  seissol::initializer::parameters::BoundaryFormat boundaryFormat;

  /**
   * Read the mesh
   */
  void read(PUML::TETPUML& puml, const char* meshFile);

  /**
   * Create the partitioning
   */
  void partition(PUML::TETPUML& puml,
                 initializer::time_stepping::LtsWeights* ltsWeights,
                 double tpwgt,
                 const char* meshFile,
                 const char* partitioningLib,
                 bool readPartitionFromFile,
                 const char* checkPointFile);
  int readPartition(PUML::TETPUML& puml, int* partition, const char* checkPointFile);
  void writePartition(PUML::TETPUML& puml, int* partition, const char* checkPointFile);
  /**
   * Generate the PUML data structure
   */
  void generatePUML(PUML::TETPUML& puml);

  /**
   * Get the mesh
   */
  void getMesh(const PUML::TETPUML& puml);

  void addMPINeighor(const PUML::TETPUML& puml, int rank, const std::vector<unsigned int>& faces);
};

} // namespace seissol::geometry

#endif // SEISSOL_SRC_GEOMETRY_PUMLREADER_H_
