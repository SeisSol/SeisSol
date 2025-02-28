// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

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
             seissol::initializer::parameters::BoundaryFormat boundaryFormat =
                 seissol::initializer::parameters::BoundaryFormat::I32,
             initializer::time_stepping::LtsWeights* ltsWeights = nullptr,
             double tpwgt = 1.0);

  bool inlineTimestepCompute() const override;
  bool inlineClusterCompute() const override;

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

#endif // SEISSOL_SRC_GEOMETRY_PUMLREADER_H_
