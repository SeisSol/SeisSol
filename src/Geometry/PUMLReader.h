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
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace seissol::initializer {
class Clustering;
struct ClusteringResult;
class VertexWeightModel;
} // namespace seissol::initializer

namespace seissol::geometry {
constexpr PUML::TopoType PumlTopology = PUML::TETRAHEDRON;
using PumlMesh = PUML::PUML<PumlTopology>;

/// the names the cell data of a mesh read by PUMLReader is filed under
namespace pumldata {
constexpr const char* Group = "group";
constexpr const char* Boundary = "boundary";
constexpr const char* CellIdsAsInFile = "cellIdsAsInFile";
constexpr const char* ClusterIds = "clusterIds";
constexpr const char* Timesteps = "timesteps";
constexpr const char* VertexOrders = "vertexOrders";
/// the order of every cell and its nodes but the vertices, where the file gives them
constexpr const char* Order = "order";
constexpr const char* HighOrderNodes = "highOrderNodes";
} // namespace pumldata

/// the group of every cell of a mesh read by PUMLReader
auto groupsOf(const PumlMesh& mesh) -> const int*;

/// the boundary condition tags of every cell of a mesh read by PUMLReader, as decodeBoundary reads
/// them
auto boundaryTagsOf(const PumlMesh& mesh, seissol::initializer::parameters::BoundaryFormat format)
    -> const void*;

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

/// The local vertex order of a cell: entry k is the slot, in the vertex list of the cell as the
/// mesh file gives it, of the vertex that SeisSol uses as local vertex k.
using VertexOrder = std::array<std::uint8_t, Cell::NumVertices>;

/**
 * The canonical local vertex order of every local cell of the two meshes, which have to
 * correspond cell by cell and slot by slot (as they do from reading until getMesh(),
 * partitioning included). It needs no communication: the sort key is the global vertex id from
 * the file.
 */
std::vector<VertexOrder> canonicalVertexOrders(const PumlMesh& meshTopology,
                                               const PumlMesh& meshGeometry);

/// What the nodes of a higher order a mesh file gives for its cells say about them, over all ranks.
struct HighOrderGeometrySummary {
  /// the cells of an order above one
  unsigned long highOrderCells{0};
  /// those of them whose nodes are not where the straight-sided cell through its vertices has them
  unsigned long curvedCells{0};
  /// the highest order of a cell
  unsigned long order{1};
  /// how far the nodes of a cell are off the straight-sided one at most, relative to its longest
  /// edge
  double curvature{0};
};

/**
 * Sums up the nodes of a higher order of the cells of a mesh read with
 * PUMLReader::readHighOrderGeometry, on all ranks; all zero where the file gives none. A cell
 * counts as curved where a node is off by more than PUMLReader::StraightTolerance of its longest
 * edge. Collective, and the mesh has to be generated.
 */
auto summarizeHighOrderGeometry(const PumlMesh& meshGeometry) -> HighOrderGeometrySummary;

class PUMLReader : public seissol::geometry::MeshReader {
  public:
  /// a cell whose nodes are this close to the straight-sided one, relative to its longest edge, is
  /// the straight-sided one for the scheme
  static constexpr double StraightTolerance = 1e-8;

  PUMLReader(const std::string& meshFile,
             const std::string& partitioningLib,
             const seissol::FaceMap& faceMap,
             seissol::initializer::parameters::BoundaryFormat boundaryFormat =
                 seissol::initializer::parameters::BoundaryFormat::I32,
             seissol::initializer::parameters::TopologyFormat topologyFormat =
                 seissol::initializer::parameters::TopologyFormat::Geometric,
             initializer::Clustering* clustering = nullptr,
             initializer::VertexWeightModel* weightModel = nullptr,
             double tpwgt = 1.0);

  bool inlineTimestepCompute() const override;
  bool inlineClusterCompute() const override;

  /**
   * Reads the order of every cell and its nodes beyond the vertices, where the file gives them, as
   * a mesh file of PUMGen does for a mesh of gmsh whose cells are of a higher order: the nodes in
   * /geometry_ho, in the order of gmsh, where those of every cell start in /geometry_ho_offsets,
   * and the order of every cell in /order. The cells have to be read already.
   */
  static void readHighOrderGeometry(PumlMesh& meshGeometry, const std::string& file);

  /**
   * Takes a mesh whose file gives nodes of a higher order only where its cells are straight-sided,
   * and then forgets the nodes, so that they do not travel with their cells through the
   * partitioning; a mesh with a curved cell it refuses, since SeisSol cannot compute on curved
   * cells yet. Collective, and the mesh has to be generated.
   */
  static void acceptHighOrderGeometry(PumlMesh& meshGeometry);

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
                        const initializer::ClusteringResult* clustering,
                        const std::vector<VertexOrder>& vertexOrders,
                        initializer::VertexWeightModel* weightModel,
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

  void addMPINeighor(const PumlMesh& meshTopology,
                     int rank,
                     const std::vector<PUML::LocalId>& faces,
                     const std::vector<std::array<std::uint8_t, Cell::NumFaces>>& pumlFaceMaps);
};

} // namespace seissol::geometry

#endif // SEISSOL_SRC_GEOMETRY_PUMLREADER_H_
