// SPDX-FileCopyrightText: 2015 SeisSol Group
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
#include "Parallel/MPI.h"

#include <PUML/PUML.h>
#include <PUML/Topology.h>

namespace seissol::initializer::time_stepping {
class LtsWeights;
} // namespace seissol::initializer::time_stepping

namespace seissol::geometry {
constexpr PUML::TopoType PumlTopology = PUML::TETRAHEDRON;
using PumlMesh = PUML::PUML<PumlTopology>;

inline uint32_t decodeBoundary(const void* data,
                               size_t cell,
                               uint8_t face,
                               seissol::initializer::parameters::BoundaryFormat format) {
  if (format == seissol::initializer::parameters::BoundaryFormat::I32) {
    const auto* dataCasted = reinterpret_cast<const uint32_t*>(data);
    return (dataCasted[cell] >> (8 * face)) & 0xff;
  } else if (format == seissol::initializer::parameters::BoundaryFormat::I64) {
    const auto* dataCasted = reinterpret_cast<const uint64_t*>(data);
    return (dataCasted[cell] >> (16 * face)) & 0xffff;
  } else if (format == seissol::initializer::parameters::BoundaryFormat::I32x4) {
    const auto* dataCasted = reinterpret_cast<const int*>(data);
    return dataCasted[cell * Cell::NumFaces + face];
  } else {
    logError() << "Unknown boundary format:" << static_cast<uint32_t>(format);
    return 0;
  }
}

inline std::optional<FaceType> boundaryTagToFaceType(uint32_t tag) {
  if (tag == 0 || tag == 6) {
    return FaceType::Regular;
  }
  if (tag == 1) {
    return FaceType::FreeSurface;
  }
  if (tag == 2) {
    return FaceType::FreeSurfaceGravity;
  }
  if (tag == 3 || tag > 64) {
    return FaceType::DynamicRupture;
  }
  if (tag == 4) {
    return FaceType::Dirichlet;
  }
  if (tag == 5) {
    return FaceType::Outflow;
  }
  if (tag == 7) {
    return FaceType::Analytical;
  }

  return {};
}

class PUMLReader : public seissol::geometry::MeshReader {
  public:
  PUMLReader(const std::string& meshFile,
             const std::string& partitioningLib,
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
               seissol::initializer::parameters::BoundaryFormat boundaryFormat);

  void
      addMPINeighor(const PumlMesh& meshTopology, int rank, const std::vector<unsigned int>& faces);
};

} // namespace seissol::geometry

#endif // SEISSOL_SRC_GEOMETRY_PUMLREADER_H_
