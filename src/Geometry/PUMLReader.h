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

#include <PUML/DataHandle.h>
#include <PUML/PUML.h>
#include <PUML/RaggedBuffer.h>
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
/// the cells by the vertices of the topology, where the file gives a topology of its own
constexpr const char* Topology = "topology";
/// the order of every cell and its nodes but the vertices, where the file gives them
constexpr const char* Order = "order";
constexpr const char* HighOrderNodes = "highOrderNodes";
} // namespace pumldata

/**
 * The arrays PUMLReader reads with the cells of a mesh, by the handles PUML hands out for them: a
 * handle knows the type of the values it names, and it reads the mesh it is of only.
 */
struct PumlCellArrays {
  PUML::DataHandle<int> group;
  /// the boundary condition tags of the faces, in `boundaryFormat` (see BoundaryTagsView)
  PUML::DataHandle<std::uint64_t> boundary;
  seissol::initializer::parameters::BoundaryFormat boundaryFormat{
      seissol::initializer::parameters::BoundaryFormat::I32};
  PUML::DataHandle<std::size_t> cellIdsAsInFile;

  /// The handles of the same arrays in `clone`, a clone of the mesh these are of.
  [[nodiscard]] auto rebasedTo(const PumlMesh& clone) const -> PumlCellArrays;
};

/**
 * The boundary condition tags of the faces of the cells, as a mesh file gives them: one value per
 * cell that packs the tags of its four faces, eight bits each in I32 and 16 in I64, or one value
 * per face in I32x4. Every value is read as a 64-bit integer whose bits are those of the file, so
 * that a tag reads the same whatever integer type the file has (see PUMLReader::read).
 */
struct BoundaryTagsView {
  PUML::DataView<const std::uint64_t> values;
  seissol::initializer::parameters::BoundaryFormat format{
      seissol::initializer::parameters::BoundaryFormat::I32};

  /// the tag of a face of a local cell, the face counted as PUML counts them
  [[nodiscard]] auto tag(std::size_t cell, std::size_t face) const -> std::uint32_t {
    using seissol::initializer::parameters::BoundaryFormat;
    if (format == BoundaryFormat::I32) {
      return static_cast<std::uint32_t>((values[cell] >> (8 * face)) & 0xff);
    }
    if (format == BoundaryFormat::I64) {
      return static_cast<std::uint32_t>((values[cell] >> (16 * face)) & 0xffff);
    }
    // I32x4, since PUMLReader::read takes no other format
    return static_cast<std::uint32_t>(values[cell * Cell::NumFaces + face]);
  }
};

/// the boundary condition tags of the faces of the cells of a mesh, which PUMLReader read as
/// `arrays`
inline auto boundaryTagsOf(const PumlMesh& mesh, const PumlCellArrays& arrays) -> BoundaryTagsView {
  return BoundaryTagsView{.values = mesh.data(arrays.boundary), .format = arrays.boundaryFormat};
}

/// The order of every cell and its nodes but the vertices, where the file gives them; neither
/// handle is valid where it gives none.
struct HighOrderArrays {
  PUML::DataHandle<std::uint8_t> order;
  PUML::RaggedHandle<double> nodes;
};

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
  std::uint64_t highOrderCells{0};
  /// those of them whose nodes are not where the straight-sided cell through its vertices has them
  std::uint64_t curvedCells{0};
  /// the highest order of a cell
  std::uint64_t order{1};
  /// how far the nodes of a cell are off the straight-sided one at most, relative to its longest
  /// edge
  double curvature{0};
};

/**
 * Sums up the nodes of a higher order of the cells of a mesh, which
 * PUMLReader::readHighOrderGeometry read as `arrays`, on all ranks; all zero where the file gives
 * none. A cell counts as curved where a node is off by more than PUMLReader::StraightTolerance of
 * its longest edge. Collective, and the mesh has to be generated.
 */
auto summarizeHighOrderGeometry(const PumlMesh& meshGeometry, const HighOrderArrays& arrays)
    -> HighOrderGeometrySummary;

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
   * Reads the cells and the vertices of a mesh file, with the arrays SeisSol keeps for the cells.
   */
  static auto read(PumlMesh& mesh,
                   const std::string& file,
                   seissol::initializer::parameters::BoundaryFormat boundaryFormat)
      -> PumlCellArrays;

  /**
   * Reads the order of every cell and its nodes beyond the vertices, where the file gives them, as
   * a mesh file of PUMGen does for a mesh of gmsh whose cells are of a higher order: the nodes in
   * /geometry_ho, in the order of gmsh, where those of every cell start in /geometry_ho_offsets,
   * and the order of every cell in /order. The cells have to be read already.
   */
  static auto readHighOrderGeometry(PumlMesh& meshGeometry, const std::string& file)
      -> HighOrderArrays;

  /**
   * Takes a mesh whose file gives nodes of a higher order only where its cells are straight-sided,
   * and then forgets the nodes, so that they do not travel with their cells through the
   * partitioning; a mesh with a curved cell it refuses, since SeisSol cannot compute on curved
   * cells yet. Collective, and the mesh has to be generated.
   */
  static void acceptHighOrderGeometry(PumlMesh& meshGeometry, const HighOrderArrays& arrays);

  private:
  /// The arrays the partitioning carries along with the cells, for getMesh.
  struct PartitionArrays {
    PUML::DataHandle<std::size_t> clusterIds;
    PUML::DataHandle<double> timesteps;
    PUML::DataHandle<std::uint8_t> vertexOrders;
  };

  /**
   * Create the partitioning
   */
  static auto partition(PumlMesh& meshTopology,
                        PumlMesh& meshGeometry,
                        const initializer::ClusteringResult* clustering,
                        const std::vector<VertexOrder>& vertexOrders,
                        initializer::VertexWeightModel* weightModel,
                        double tpwgt,
                        const std::string& partitioningLib) -> PartitionArrays;
  /**
   * Generate the PUML data structure
   */
  static void generatePUML(PumlMesh& meshTopology, PumlMesh& meshGeometry);

  /**
   * Get the mesh
   */
  void getMesh(const PumlMesh& meshTopology,
               const PumlMesh& meshGeometry,
               const PumlCellArrays& cellArrays,
               const PartitionArrays& partitionArrays,
               const FaceMap& faceMap);

  void addMPINeighor(const PumlMesh& meshTopology,
                     int rank,
                     const std::vector<PUML::LocalId>& faces,
                     const std::vector<std::array<std::uint8_t, Cell::NumFaces>>& pumlFaceMaps);
};

} // namespace seissol::geometry

#endif // SEISSOL_SRC_GEOMETRY_PUMLREADER_H_
