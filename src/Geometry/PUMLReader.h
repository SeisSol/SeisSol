// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_GEOMETRY_PUMLREADER_H_
#define SEISSOL_SRC_GEOMETRY_PUMLREADER_H_

#include "Initializer/FaceMap.h"
#include "Initializer/Parameters/MeshParameters.h"
#include "MeshReader.h"
#include "Parallel/MPI.h"

#include <PUML/PUML.h>
#include <PUML/Topology.h>

namespace seissol::initializer::time_stepping {
class LtsWeights;
} // namespace seissol::initializer::time_stepping

namespace seissol::geometry {
constexpr PUML::TopoType PumlTopology = PUML::TETRAHEDRON;
using PumlMesh = PUML::PUML<PumlTopology>;

/**
 * The boundary condition tag of a face, from the tags as PUMLReader reads them: as 64-bit integers
 * with the bits they have in the file, in all formats.
 */
inline uint32_t decodeBoundary(const void* data,
                               size_t cell,
                               uint8_t face,
                               seissol::initializer::parameters::BoundaryFormat format) {
  const auto* dataCasted = reinterpret_cast<const uint64_t*>(data);
  if (format == seissol::initializer::parameters::BoundaryFormat::I32) {
    return (dataCasted[cell] >> (8 * face)) & 0xff;
  } else if (format == seissol::initializer::parameters::BoundaryFormat::I64) {
    return (dataCasted[cell] >> (16 * face)) & 0xffff;
  } else if (format == seissol::initializer::parameters::BoundaryFormat::I32x4) {
    return static_cast<uint32_t>(dataCasted[cell * Cell::NumFaces + face]);
  } else {
    logError() << "Unknown boundary format:" << static_cast<uint32_t>(format);
    return 0;
  }
}

class PUMLReader : public seissol::geometry::MeshReader {
  public:
  PUMLReader(const std::string& meshFile,
             const std::string& partitioningLib,
             const seissol::FaceMap& faceMap,
             seissol::initializer::parameters::BoundaryFormat boundaryFormat =
                 seissol::initializer::parameters::BoundaryFormat::I32,
             seissol::initializer::parameters::TopologyFormat topologyFormat =
                 seissol::initializer::parameters::TopologyFormat::Geometric,
             initializer::time_stepping::LtsWeights* ltsWeights = nullptr,
             double tpwgt = 1.0);

  bool inlineTimestepCompute() const override;
  bool inlineClusterCompute() const override;

  private:
  /**
   * Read the mesh
   */
  static void read(PumlMesh& meshTopology,
                   const std::string& file,
                   bool topology,
                   seissol::initializer::parameters::BoundaryFormat boundaryFormat);

  /**
   * Create the partitioning
   */
  static void partition(PumlMesh& meshTopology,
                        PumlMesh& meshGeometry,
                        initializer::time_stepping::LtsWeights* ltsWeights,
                        double tpwgt,
                        const std::string& partitioningLib);
  /**
   * Generate the PUML data structure
   */
  static void generatePUML(PumlMesh& meshTopology, PumlMesh& meshGeometry);

  /**
   * Get the mesh
   */
  void getMesh(const PumlMesh& meshTopology,
               const PumlMesh& meshGeometry,
               const FaceMap& faceMap,
               seissol::initializer::parameters::BoundaryFormat boundaryFormat);

  void
      addMPINeighor(const PumlMesh& meshTopology, int rank, const std::vector<unsigned int>& faces);
};

} // namespace seissol::geometry

#endif // SEISSOL_SRC_GEOMETRY_PUMLREADER_H_
