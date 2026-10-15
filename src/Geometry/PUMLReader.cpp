// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#include "PUMLReader.h"

#include "Common/CompactOptional.h"
#include "Common/Constants.h"
#include "Common/Iterator.h"
#include "Geometry/CellTransform.h"
#include "Geometry/GmshNodes.h"
#include "Geometry/MeshDefinition.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/Clustering/Clustering.h"
#include "Initializer/Clustering/VertexWeights/VertexWeightModel.h"
#include "Initializer/FaceMap.h"
#include "Initializer/Parameters/MeshParameters.h"
#include "Monitoring/Instrumentation.h"
#include "PartitioningLib.h"

#include <Eigen/Dense>
#include <PUML/DataHandle.h>
#include <PUML/Downward.h>
#include <PUML/Hdf5Reader.h>
#include <PUML/Neighbor.h>
#include <PUML/PUML.h>
#include <PUML/Partition.h>
#include <PUML/PartitionGraph.h>
#include <PUML/PartitionTarget.h>
#include <PUML/RaggedBuffer.h>
#include <PUML/Topology.h>
#include <PUML/TypeInference.h>
#include <PUML/Types.h>
#include <PUML/Upward.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <hdf5.h>
#include <mpi.h>
#include <numeric>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utils/logger.h>
#include <vector>

namespace seissol::geometry {

namespace {

// PUML sanity checks
using PumlTopologyInternal = PUML::internal::Topology<seissol::geometry::PumlTopology>;
static_assert(PumlTopologyInternal::dimension() == seissol::Cell::Dim);
static_assert(PumlTopologyInternal::cellfaces() == seissol::Cell::NumFaces);
static_assert(PumlTopologyInternal::cellvertices() == seissol::Cell::NumVertices);

void logassertI(bool condition, const std::string& file, int line) {
  if (!condition) {
    logError() << "Assertion failure in" << file << "at" << line;
  }
}

#define logassert(x) logassertI(x, __FILE__, __LINE__)

/**
 * Decodes the boundary condition tag into a string representation.
 */
inline std::string bcToString(uint32_t id, const FaceMap& faceMap) {
  const auto type = faceMap.at(id);
  if (type == FaceType::Regular) {
    return std::string("regular");
  } else if (type == FaceType::FreeSurface) {
    return std::string("free surface");
  } else if (type == FaceType::FreeSurfaceGravity) {
    return std::string("free surface with gravity");
  } else if (type == FaceType::Dirichlet) {
    return std::string("dirichlet");
  } else if (type == FaceType::Outflow) {
    return std::string("outflow");
  } else if (type == FaceType::Analytical) {
    return std::string("analytical");
  } else if (type == FaceType::DynamicRupture) {
    std::stringstream s;
    s << "dynamic rupture (face tag " << id << ")";
    return s.str();
  } else {
    std::stringstream s;
    s << "unknown (" << id << ")";
    return s.str();
  }
}

/**
 * Check if the mesh is locally correct:
 * - if a face is an internal face, the neighbor has to exist.
 * - if a face is an external face, the neighbor must not exist.
 * @param face: face to check
 * @param cellNeighbors: ids of the neighboring cells
 * @param side: side of the tetrahedron to check
 * @param sideBC: boundary condition tag at the side to check
 * @param cellIdAsInFile: Original cell id as it is given in the h5 file
 */
template <PUML::TopoType Topo>
inline bool checkMeshCorrectnessLocally(
    const typename PUML::PUML<Topo>::face_t& face,
    const std::array<PUML::LocalId, seissol::Cell::NumFaces>& cellNeighbors,
    uint8_t side,
    uint32_t sideBC,
    uint64_t cellIdAsInFile,
    const FaceMap& faceMap) {
  // all of these will only issue warnings here -- the "logError()" is supposed to come later, after
  // all warning have been logged

  const auto faceType = faceMap.at(sideBC);

  if (faceType.has_value()) {

    // if a face is an internal face, it has to have a neighbor on either this rank or somewhere
    // else:
    if (getBCType(faceType.value()) == BCType::Internal) {
      if (cellNeighbors[side] == PUML::InvalidLocalId && !face.isShared()) {
        logWarning(true) << "Element" << cellIdAsInFile << ", side" << side << " has a"
                         << bcToString(sideBC, faceMap)
                         << "boundary condition, but the neighboring element doesn't exist";
        return false;
      }
    }
    // external boundaries must not have neighboring elements:
    else {
      if (cellNeighbors[side] != PUML::InvalidLocalId || face.isShared()) {
        logWarning(true) << "Element" << cellIdAsInFile << ", side" << side << " has a"
                         << bcToString(sideBC, faceMap)
                         << "boundary condition, but a neighboring element exists";
        return false;
      }
    }
  } else {
    // ignore unknown boundary conditions and warn
    logWarning(true) << "Element" << cellIdAsInFile << ", side" << side
                     << " has a boundary condition (" << sideBC
                     << ") which is not understood by this version of SeisSol";
    return false;
  }
  return true;
}

// helper arrays

// indexes the vertices on each face i (or FaceVertexToOrientation[i][j] == -1 to indicate that the
// vertex does not lie on it)
const std::array<std::array<std::int32_t, 4>, 4> FaceVertexToOrientation = {
    std::array<std::int32_t, 4>{0, 2, 1, -1},
    std::array<std::int32_t, 4>{0, 1, -1, 2},
    std::array<std::int32_t, 4>{0, -1, 2, 1},
    std::array<std::int32_t, 4>{-1, 0, 1, 2}};

// the first vertex on the face (i.e. FirstFaceVertex[i] == j, where j is the lowest index in
// FaceVertexToOrientation[i] to not be -1)
const std::array<std::int32_t, 4> FirstFaceVertex = {0, 0, 0, 1};

// PUML face p is opposite the local vertex PumlFaceMissingVertex[p]; SeisSol face s is opposite
// the local vertex Cell::NumVertices - 1 - s
const std::array<std::size_t, 4> PumlFaceMissingVertex = {3, 2, 0, 1};

/**
 * Signed volume measure of a tetrahedron, with the same convention as
 * MeshReader::verifyMeshOrientation, which requires a negative value.
 */
double orientationDeterminant(const std::array<const double*, Cell::NumVertices>& coords) {
  Eigen::Matrix<double, 4, 4> mat;
  mat << coords[0][0], coords[0][1], coords[0][2], 1, coords[1][0], coords[1][1], coords[1][2], 1,
      coords[2][0], coords[2][1], coords[2][2], 1, coords[3][0], coords[3][1], coords[3][2], 1;
  return mat.determinant();
}

/**
 * The canonical local vertex numbering of a cell: sort the vertices by global topological id,
 * then repair the orientation with the transposition (2 3) where the sorted order would come out
 * positively oriented. Returns the map from the new local slot to the old one.
 *
 * Sorting alone already forces the face orientation index to zero on every interior face.
 * FirstFaceVertex is the smallest local index on each face, hence under a sorted numbering it is
 * the vertex with the smallest id on that face -- and both sides of a face agree on which vertex
 * that is. The transposition (2 3) leaves local slots 0 and 1 alone and therefore keeps that
 * property intact, which no other odd permutation does.
 *
 * The sort key has to be the global topological id: the topology identifies the two sides of a
 * face, and only the global id is an ordering that both sides agree on across rank boundaries,
 * where the rank-local vertex indices are unrelated to each other. The orientation in turn has to
 * be judged on the geometry. For periodic meshes the two differ.
 */
VertexOrder canonicalVertexOrder(const std::array<unsigned long, Cell::NumVertices>& topoVertices,
                                 const std::array<const double*, Cell::NumVertices>& coords) {
  VertexOrder order{};
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](auto a, auto b) {
    return topoVertices[a] < topoVertices[b];
  });

  std::array<const double*, Cell::NumVertices> sorted{};
  for (std::size_t k = 0; k < Cell::NumVertices; ++k) {
    sorted[k] = coords[order[k]];
  }
  if (orientationDeterminant(sorted) >= 0) {
    std::swap(order[2], order[3]);
  }
  return order;
}

/**
 * PUML face index -> SeisSol face index for a given canonical vertex order. Both index schemes
 * cover the same four faces; only the ordering within a cell differs.
 */
std::array<std::uint8_t, Cell::NumFaces> canonicalFaceMap(const VertexOrder& order) {
  std::array<std::size_t, Cell::NumVertices> inverse{};
  for (std::size_t k = 0; k < Cell::NumVertices; ++k) {
    inverse[order[k]] = k;
  }
  std::array<std::uint8_t, Cell::NumFaces> map{};
  for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
    map[f] = static_cast<std::uint8_t>(Cell::NumVertices - 1 - inverse[PumlFaceMissingVertex[f]]);
  }
  return map;
}

/**
 * Recover the vertex order from a face map. PumlFaceMissingVertex is a bijection, so the face map
 * already determines the permutation and only one of the two needs to be kept per cell.
 */
VertexOrder vertexOrderFromFaceMap(const std::array<std::uint8_t, Cell::NumFaces>& map) {
  VertexOrder order{};
  for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
    order[Cell::NumVertices - 1 - map[f]] = static_cast<std::uint8_t>(PumlFaceMissingVertex[f]);
  }
  return order;
}

template <typename T>
void applyVertexOrder(std::array<T, Cell::NumVertices>& values, const VertexOrder& order) {
  const auto original = values;
  for (std::size_t k = 0; k < Cell::NumVertices; ++k) {
    values[k] = original[order[k]];
  }
}

/**
 * A string attribute of a dataset, read on the first rank and handed to the others; empty where
 * the dataset has no such attribute.
 */
std::string stringAttribute(const std::string& file, const char* dataset, const char* name) {
  std::string value;
  if (Mpi::mpi.rank() == 0) {
    const hid_t h5file = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (h5file >= 0) {
      if (H5Aexists_by_name(h5file, dataset, name, H5P_DEFAULT) > 0) {
        const hid_t attribute = H5Aopen_by_name(h5file, dataset, name, H5P_DEFAULT, H5P_DEFAULT);
        const hid_t fileType = H5Aget_type(attribute);
        if (H5Tget_class(fileType) == H5T_STRING) {
          const hid_t memType = H5Tcopy(H5T_C_S1);
          if (H5Tis_variable_str(fileType) > 0) {
            H5Tset_size(memType, H5T_VARIABLE);
            char* text = nullptr;
            if (H5Aread(attribute, memType, static_cast<void*>(&text)) >= 0 && text != nullptr) {
              value = text;
              H5free_memory(text);
            }
          } else {
            std::vector<char> text(H5Tget_size(fileType) + 1, '\0');
            H5Tset_size(memType, text.size());
            if (H5Aread(attribute, memType, text.data()) >= 0) {
              value = text.data();
            }
          }
          H5Tclose(memType);
        }
        H5Tclose(fileType);
        H5Aclose(attribute);
      }
      H5Fclose(h5file);
    }
  }
  auto length = static_cast<unsigned long>(value.size());
  MPI_Bcast(&length, 1, MPI_UNSIGNED_LONG, 0, Mpi::mpi.comm());
  value.resize(length);
  MPI_Bcast(value.data(), static_cast<int>(length), MPI_CHAR, 0, Mpi::mpi.comm());
  return value;
}

/**
 * Whether a dataset holds signed integers (1), unsigned ones (0), or anything else, integers wider
 * than 64 bits included (-1), as the first rank finds; all ranks learn it.
 */
int integerSignOf(const std::string& file, const char* dataset) {
  int sign = -1;
  if (Mpi::mpi.rank() == 0) {
    const hid_t h5file = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (h5file >= 0) {
      const hid_t h5dataset = H5Dopen2(h5file, dataset, H5P_DEFAULT);
      if (h5dataset >= 0) {
        const hid_t fileType = H5Dget_type(h5dataset);
        if (H5Tget_class(fileType) == H5T_INTEGER &&
            H5Tget_size(fileType) <= sizeof(std::uint64_t)) {
          sign = H5Tget_sign(fileType) == H5T_SGN_2 ? 1 : 0;
        }
        H5Tclose(fileType);
        H5Dclose(h5dataset);
      }
      H5Fclose(h5file);
    }
  }
  MPI_Bcast(&sign, 1, MPI_INT, 0, Mpi::mpi.comm());
  return sign;
}

/**
 * Where the cells a rank reads start among the cells of the file, which the ranks read in runs, one
 * after the other: the number of cells the ranks before it read. Collective.
 */
std::size_t firstCellOf(const PumlMesh& mesh) {
  const std::size_t localCells = mesh.numOriginalCells();
  std::size_t firstCell = 0;
  MPI_Exscan(
      &localCells, &firstCell, 1, PUML::MPITypeInfer<std::size_t>::type(), MPI_SUM, mesh.comm());
  // MPI_Exscan leaves the result on the first rank undefined
  int rank = 0;
  MPI_Comm_rank(mesh.comm(), &rank);
  return rank == 0 ? 0 : firstCell;
}

/// the vertices of a local cell, in the vertex order of the file
std::array<CellTransform::VectorEigenT, Cell::NumVertices>
    fileVertices(const PumlMesh& meshGeometry, std::size_t cell) {
  std::array<PUML::LocalId, Cell::NumVertices> vertexIds{};
  PUML::Downward::vertices(meshGeometry, meshGeometry.cells()[cell], vertexIds.data());
  std::array<CellTransform::VectorEigenT, Cell::NumVertices> vertices;
  for (std::size_t k = 0; k < Cell::NumVertices; ++k) {
    const auto* coordinate = meshGeometry.vertices()[vertexIds[k]].coordinate();
    vertices[k] = CellTransform::VectorEigenT(coordinate[0], coordinate[1], coordinate[2]);
  }
  return vertices;
}

/// the order of every local cell and its nodes but the vertices, where the file gives them
class HighOrderGeometry {
  public:
  HighOrderGeometry(const PumlMesh& meshGeometry, const HighOrderArrays& arrays)
      : present_(arrays.order.valid()) {
    if (present_) {
      orders_ = meshGeometry.data(arrays.order);
      nodes_ = meshGeometry.raggedData(arrays.nodes);
    }
  }

  [[nodiscard]] bool present() const { return present_; }

  [[nodiscard]] std::size_t order(std::size_t cell) const { return present_ ? orders_[cell] : 1; }

  /// the nodes of a cell but its vertices, in the order of gmsh, three coordinates each
  [[nodiscard]] const double* nodes(std::size_t cell) const {
    return present_ ? nodes_.begin(cell) : nullptr;
  }

  /// how many values the file gives for a cell
  [[nodiscard]] std::size_t values(std::size_t cell) const {
    return present_ ? nodes_.count(cell) : 0;
  }

  private:
  bool present_;
  PUML::DataView<const std::uint8_t> orders_;
  PUML::RaggedView<const double> nodes_;
};

/**
 * How far the nodes of a local cell are off the straight-sided cell through its vertices, relative
 * to its longest edge: the straight-sided cell has the node gmsh numbers i where the affine map of
 * the cell puts the point of the reference cell that node stands for.
 */
double
    curvatureOf(const PumlMesh& meshGeometry, const HighOrderGeometry& geometry, std::size_t cell) {
  const auto order = geometry.order(cell);
  if (order <= 1) {
    return 0;
  }
  const auto vertices = fileVertices(meshGeometry, cell);
  const AffineTransform straight(vertices);
  double size = 0;
  for (std::size_t a = 0; a < Cell::NumVertices; ++a) {
    for (std::size_t b = a + 1; b < Cell::NumVertices; ++b) {
      size = std::max(size, (vertices[a] - vertices[b]).norm());
    }
  }
  const auto lattice = gmshTetrahedronLattice(order);
  const auto* nodes = geometry.nodes(cell);
  double offset = 0;
  for (std::size_t i = Cell::NumVertices; i < lattice.size(); ++i) {
    const CellTransform::VectorEigenT reference(static_cast<double>(lattice[i][0]),
                                                static_cast<double>(lattice[i][1]),
                                                static_cast<double>(lattice[i][2]));
    const auto* node = nodes + Cell::Dim * (i - Cell::NumVertices);
    const CellTransform::VectorEigenT given(node[0], node[1], node[2]);
    offset = std::max(offset,
                      (given - straight.refToSpace(reference / static_cast<double>(order))).norm());
  }
  return offset / size;
}
} // namespace

auto PumlCellArrays::rebasedTo(const PumlMesh& clone) const -> PumlCellArrays {
  auto rebased = *this;
  rebased.group = clone.rebase(group);
  rebased.boundary = clone.rebase(boundary);
  rebased.cellIdsAsInFile = clone.rebase(cellIdsAsInFile);
  return rebased;
}

auto summarizeHighOrderGeometry(const PumlMesh& meshGeometry, const HighOrderArrays& arrays)
    -> HighOrderGeometrySummary {
  const HighOrderGeometry geometry(meshGeometry, arrays);
  HighOrderGeometrySummary summary;
  if (geometry.present()) {
    for (std::size_t cell = 0; cell < meshGeometry.cells().size(); ++cell) {
      if (geometry.order(cell) > 1) {
        ++summary.highOrderCells;
        summary.order = std::max<unsigned long>(summary.order, geometry.order(cell));
        const auto curvature = curvatureOf(meshGeometry, geometry, cell);
        summary.curvature = std::max(summary.curvature, curvature);
        summary.curvedCells += curvature > PUMLReader::StraightTolerance ? 1 : 0;
      }
    }
  }
  MPI_Allreduce(
      MPI_IN_PLACE, &summary.highOrderCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, Mpi::mpi.comm());
  MPI_Allreduce(MPI_IN_PLACE, &summary.curvedCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, Mpi::mpi.comm());
  MPI_Allreduce(MPI_IN_PLACE, &summary.order, 1, MPI_UNSIGNED_LONG, MPI_MAX, Mpi::mpi.comm());
  MPI_Allreduce(MPI_IN_PLACE, &summary.curvature, 1, MPI_DOUBLE, MPI_MAX, Mpi::mpi.comm());
  return summary;
}

std::vector<VertexOrder> canonicalVertexOrders(const PumlMesh& meshTopology,
                                               const PumlMesh& meshGeometry) {
  const std::vector<PumlMesh::cell_t>& cells = meshTopology.cells();
  const std::vector<PumlMesh::cell_t>& cellsGeometry = meshGeometry.cells();
  const std::vector<PumlMesh::vertex_t>& verticesGeometry = meshGeometry.vertices();
  assert(cells.size() == cellsGeometry.size());

  std::vector<VertexOrder> orders(cells.size());
  for (std::size_t i = 0; i < cells.size(); i++) {
    std::array<unsigned long, Cell::NumVertices> topoVertices{};
    PUML::Downward::gvertices(meshTopology, cells[i], topoVertices.data());

    std::array<PUML::LocalId, Cell::NumVertices> geomVertices{};
    PUML::Downward::vertices(meshGeometry, cellsGeometry[i], geomVertices.data());

    std::array<const double*, Cell::NumVertices> coords{};
    for (std::size_t k = 0; k < Cell::NumVertices; k++) {
      coords[k] = verticesGeometry[geomVertices[k]].coordinate();
    }

    orders[i] = canonicalVertexOrder(topoVertices, coords);
  }
  return orders;
}

PUMLReader::PUMLReader(const std::string& meshFile,
                       const std::string& partitioningLib,
                       const seissol::FaceMap& faceMap,
                       seissol::initializer::parameters::BoundaryFormat boundaryFormat,
                       seissol::initializer::parameters::TopologyFormat topologyFormat,
                       initializer::Clustering* clustering,
                       initializer::VertexWeightModel* weightModel,
                       double tpwgt) {
  // we need up to two meshes, potentially:
  // one mesh for the geometry
  // one mesh for the topology
  // they will only differ if we have periodic boundary conditions

  PumlMesh meshTopologyExtra;
  PumlMesh meshGeometry;
  meshTopologyExtra.setComm(seissol::Mpi::mpi.comm());
  meshGeometry.setComm(seissol::Mpi::mpi.comm());

  logInfo() << "Read (geometric) connectivity data.";
  const auto cellArrays = read(meshGeometry, meshFile, boundaryFormat);
  const auto highOrderArrays = readHighOrderGeometry(meshGeometry, meshFile);

  // Note: we need to call generatePUML in order to create the dual graph of the mesh
  // Note 2: we also need it for vertex identification
  logInfo() << "Generate (geometric) mesh.";
  meshGeometry.generateMesh();
  acceptHighOrderGeometry(meshGeometry, highOrderArrays);

  // the arrays of the cells in the topology mesh
  auto topologyArrays = cellArrays;
  if (topologyFormat != initializer::parameters::TopologyFormat::Geometric) {
    // We have a topology mesh, separate from the physical one. It holds the same cells, with the
    // same data, and the same vertices, so it is a clone of the physical mesh rather than read from
    // the file once more; only what makes it a topology of its own comes from there.
    logInfo() << "Clone the (geometric) mesh for the topology.";
    meshTopologyExtra = meshGeometry.clone();
    topologyArrays = cellArrays.rebasedTo(meshTopologyExtra);

    PUML::DataHandle<PUML::GlobalId> identify;
    {
      PUML::Hdf5Reader<PumlTopology> reader(meshTopologyExtra);
      if (topologyFormat == initializer::parameters::TopologyFormat::IdentifyFace) {
        // the cells by the vertices of the topology, in place of those of the geometry
        logInfo() << "Read (topologic) connectivity data.";
        meshTopologyExtra.setConnectivity(reader.addData<PUML::GlobalId>(pumldata::Topology,
                                                                         meshFile + ":/topology",
                                                                         PUML::DataType::Cell,
                                                                         {Cell::NumVertices}));
      }
      if (topologyFormat == initializer::parameters::TopologyFormat::IdentifyVertex) {
        logInfo() << "Read topologic identification data.";
        identify = reader.addData<PUML::GlobalId>(
            "identify", meshFile + ":/identify", PUML::DataType::Vertex, {});
      }
    }

    // generate the topology mesh for the dual graph
    logInfo() << "Generate (topologic) mesh.";
    meshTopologyExtra.generateMesh();

    if (topologyFormat == initializer::parameters::TopologyFormat::IdentifyVertex) {
      // re-identify vertices; then re-distribute
      logInfo() << "Identify topologic vertex data.";
      meshTopologyExtra.identify(meshTopologyExtra.connectivity(), identify);

      logInfo() << "Generate (topologic) mesh again.";
      meshTopologyExtra.generateMesh();
    }
  }

  auto& meshTopology = topologyFormat == initializer::parameters::TopologyFormat::Geometric
                           ? meshGeometry
                           : meshTopologyExtra;

  // Everything that looks at the vertices of a cell from here on -- the clustering first, getMesh
  // later -- has to see them in the canonical order. It is decided once, here, and carried through
  // the partitioning with the cell.
  const auto vertexOrders = canonicalVertexOrders(meshTopology, meshGeometry);

  // The clustering needs the meshes, which only exist here -- hence the orchestrator is passed
  // in and run rather than its result.
  const initializer::ClusteringResult* clusteringResult = nullptr;
  if (clustering != nullptr) {
    logInfo() << "Compute clustering.";
    clusteringResult =
        &clustering->compute(meshTopology, topologyArrays, meshGeometry, cellArrays, vertexOrders);
  }

  logInfo() << "Partition the mesh.";
  const auto partitionArrays = partition(meshTopology,
                                         meshGeometry,
                                         clusteringResult,
                                         vertexOrders,
                                         weightModel,
                                         tpwgt,
                                         partitioningLib);

  logInfo() << "Generate the correctly-distributed meshes.";
  generatePUML(meshTopology, meshGeometry);

  logInfo() << "Set up mesh data structures.";
  getMesh(meshTopology, meshGeometry, cellArrays, partitionArrays, faceMap);
}

auto PUMLReader::read(PumlMesh& mesh,
                      const std::string& file,
                      seissol::initializer::parameters::BoundaryFormat boundaryFormat)
    -> PumlCellArrays {
  SCOREP_USER_REGION("PUMLReader_read", SCOREP_USER_REGION_TYPE_FUNCTION);

  using seissol::initializer::parameters::BoundaryFormat;
  if (boundaryFormat != BoundaryFormat::I32 && boundaryFormat != BoundaryFormat::I64 &&
      boundaryFormat != BoundaryFormat::I32x4) {
    logError() << "Unknown boundary format:" << static_cast<int>(boundaryFormat);
  }

  // The boundary condition tags are read as 64-bit integers, which hold every format. HDF5
  // converts each value to the type it is read as, clipped to the range of that type, and PUMGen
  // files the tags as signed integers: in I32 and I64, a tag of 128 resp. 32768 or more on the last
  // face of a cell makes the value of the cell negative, which an unsigned type would read as 0 --
  // all four tags of the cell lost. As a 64-bit integer of the signedness of the file, a value
  // comes out sign- or zero-extended, with the bits it has in the file.
  const auto boundarySign = integerSignOf(file, "/boundary");
  if (boundarySign < 0) {
    logError() << "The boundary condition tags in" << file
               << "(/boundary) are no integers of 64 bits or less.";
  }

  PumlCellArrays arrays;
  arrays.boundaryFormat = boundaryFormat;

  PUML::Hdf5Reader<PumlTopology> reader(mesh);
  reader.open(file + ":/connect", file + ":/geometry");
  arrays.group = reader.addData<int>(pumldata::Group, file + ":/group", PUML::DataType::Cell, {});
  arrays.boundary = reader.addData<std::uint64_t>(
      pumldata::Boundary,
      file + ":/boundary",
      PUML::DataType::Cell,
      boundaryFormat == BoundaryFormat::I32x4 ? std::vector<std::size_t>{Cell::NumFaces}
                                              : std::vector<std::size_t>{},
      PUML::MPITypeInfer<std::uint64_t>::type(),
      boundarySign > 0 ? H5T_NATIVE_INT64 : H5T_NATIVE_UINT64);

  std::vector<std::size_t> cellIdsAsInFile(mesh.numOriginalCells());
  std::iota(cellIdsAsInFile.begin(), cellIdsAsInFile.end(), firstCellOf(mesh));
  arrays.cellIdsAsInFile = mesh.addDataArray(
      pumldata::CellIdsAsInFile, cellIdsAsInFile.data(), PUML::DataType::Cell, {});
  return arrays;
}

auto PUMLReader::readHighOrderGeometry(PumlMesh& meshGeometry, const std::string& file)
    -> HighOrderArrays {
  HighOrderArrays arrays;
  {
    PUML::Hdf5Reader<PumlTopology> reader(meshGeometry);
    if (!reader.exists(file + ":/geometry_ho")) {
      return arrays;
    }
    logInfo() << "Read the nodes of the cells of a higher order.";
    arrays.order =
        reader.addData<std::uint8_t>(pumldata::Order, file + ":/order", PUML::DataType::Cell, {});
    arrays.nodes = reader.readRaggedData<double>(
        pumldata::HighOrderNodes, file + ":/geometry_ho", file + ":/geometry_ho_offsets");
  }

  const auto ordering = stringAttribute(file, "/geometry_ho", "node-ordering");
  if (ordering != "gmsh") {
    logError() << "The nodes of the cells of a higher order are given in the order"
               << (ordering.empty() ? std::string("(none given)") : ordering)
               << "-- the reader knows the one of gmsh only";
  }

  const std::size_t localCells = meshGeometry.numOriginalCells();
  const auto firstCell = firstCellOf(meshGeometry);

  const HighOrderGeometry geometry(meshGeometry, arrays);
  for (std::size_t cell = 0; cell < localCells; ++cell) {
    const auto order = geometry.order(cell);
    if (order < 1) {
      logError() << "The cell" << firstCell + cell << "has the order" << order;
    }
    const auto expected = Cell::Dim * (lagrangeNodeCount(order) - Cell::NumVertices);
    if (geometry.values(cell) != expected) {
      logError() << "The cell" << firstCell + cell << "of order" << order << "has"
                 << geometry.values(cell)
                 << "coordinates of its nodes beyond the vertices, where it needs" << expected;
    }
  }
  return arrays;
}

void PUMLReader::acceptHighOrderGeometry(PumlMesh& meshGeometry, const HighOrderArrays& arrays) {
  if (!arrays.order.valid()) {
    return;
  }
  const auto summary = summarizeHighOrderGeometry(meshGeometry, arrays);
  logInfo() << "The mesh gives nodes of an order up to" << summary.order << "for"
            << summary.highOrderCells << "cells, of which" << summary.curvedCells
            << "are curved, their nodes off the straight-sided cell by up to" << summary.curvature
            << "of its longest edge.";
  if (summary.curvedCells > 0) {
    logError() << "The mesh has" << summary.curvedCells
               << "curved cells, which SeisSol cannot compute on yet. Mesh it with cells of the "
                  "first order, or with straight-sided cells.";
  }
  // the cells are straight-sided, and their nodes need not travel with them
  meshGeometry.allocateRaggedData<double>(
      pumldata::HighOrderNodes, std::vector<PUML::Size>(meshGeometry.numOriginalCells(), 0));
}

auto PUMLReader::partition(PumlMesh& meshTopology,
                           PumlMesh& meshGeometry,
                           const initializer::ClusteringResult* clustering,
                           const std::vector<VertexOrder>& vertexOrders,
                           initializer::VertexWeightModel* weightModel,
                           double tpwgt,
                           const std::string& partitioningLib) -> PartitionArrays {
  SCOREP_USER_REGION("PUMLReader_partition", SCOREP_USER_REGION_TYPE_FUNCTION);

  auto partType = toPartitionerType(std::string_view(partitioningLib));
  logInfo() << "Using the" << toStringView(partType) << "partition library and strategy.";
  if (partType == PUML::PartitionerType::None && Mpi::mpi.size() > 1) {
    logWarning()
        << partitioningLib
        << "not found. The performance of this run will probably be much lower than it could be.";
  }
  auto partitioner = PUML::TETPartition::getPartitioner(partType);
  if (partitioner == nullptr) {
    logError() << "Unrecognized partition library: " << partitioningLib;
  }
  auto graph = PUML::TETPartitionGraph(meshTopology);
  weightModel->build(*clustering);
  graph.setVertexWeights(weightModel->vertexWeights(), weightModel->nWeightsPerVertex());

  auto nodeWeights = std::vector<double>(Mpi::mpi.size());
  MPI_Allgather(&tpwgt, 1, MPI_DOUBLE, nodeWeights.data(), 1, MPI_DOUBLE, seissol::Mpi::mpi.comm());
  double sum = 0.0;
  for (const auto& w : nodeWeights) {
    sum += w;
  }
  for (auto& w : nodeWeights) {
    w /= sum;
  }

  auto target = PUML::PartitionTarget{};
  target.setPartitionWeights(nodeWeights);
  target.setImbalance(weightModel->imbalances()[0] - 1.0);

  auto newPartition = partitioner->partition(graph, target);

  // carried along with the cells, for getMesh
  PartitionArrays arrays;
  arrays.clusterIds = meshGeometry.addDataArray(
      pumldata::ClusterIds, clustering->clusterIds.data(), PUML::DataType::Cell, {});
  arrays.timesteps = meshGeometry.addDataArray(pumldata::Timesteps,
                                               clustering->timesteps.cellTimeStepWidths.data(),
                                               PUML::DataType::Cell,
                                               {});

  // the vertex order the clustering has used
  static_assert(sizeof(VertexOrder) == Cell::NumVertices * sizeof(std::uint8_t));
  arrays.vertexOrders =
      meshGeometry.addDataArray(pumldata::VertexOrders,
                                reinterpret_cast<const std::uint8_t*>(vertexOrders.data()),
                                PUML::DataType::Cell,
                                {Cell::NumVertices});

  meshGeometry.partition(newPartition.data());
  if (&meshTopology != &meshGeometry) {
    meshTopology.partition(newPartition.data());
  }
  return arrays;
}

void PUMLReader::generatePUML(PumlMesh& meshTopology, PumlMesh& meshGeometry) {
  SCOREP_USER_REGION("PUMLReader_generate", SCOREP_USER_REGION_TYPE_FUNCTION);

  if (&meshTopology != &meshGeometry) {
    logInfo() << "Generate the correct (topologic) mesh.";
    meshTopology.generateMesh();
  }
  logInfo() << "Generate the correct (geometric) mesh.";
  meshGeometry.generateMesh();
}

void PUMLReader::getMesh(const PumlMesh& meshTopology,
                         const PumlMesh& meshGeometry,
                         const PumlCellArrays& cellArrays,
                         const PartitionArrays& partitionArrays,
                         const FaceMap& faceMap) {
  SCOREP_USER_REGION("PUMLReader_getmesh", SCOREP_USER_REGION_TYPE_FUNCTION);

  const int rank = Mpi::mpi.rank();

  const std::vector<PumlMesh::cell_t>& cells = meshTopology.cells();
  const std::vector<PumlMesh::face_t>& faces = meshTopology.faces();
  const std::vector<PumlMesh::vertex_t>& vertices = meshTopology.vertices();

  const std::vector<PumlMesh::cell_t>& cellsGeometry = meshGeometry.cells();
  const std::vector<PumlMesh::vertex_t>& verticesGeometry = meshGeometry.vertices();

  const auto group = meshGeometry.data(cellArrays.group);
  const auto boundaryTags = boundaryTagsOf(meshGeometry, cellArrays);
  const auto cellIdsAsInFile = meshGeometry.data(cellArrays.cellIdsAsInFile);
  const auto clusterIds = meshGeometry.data(partitionArrays.clusterIds);
  const auto timestep = meshGeometry.data(partitionArrays.timesteps);
  const auto vertexOrders = meshGeometry.data(partitionArrays.vertexOrders);

  // List of shared local face ids
  std::unordered_map<int, std::vector<PUML::LocalId>> neighborInfo;

  bool isMeshCorrect = true;

  // Canonical local vertex numbering, as decided before the partitioning. Mapped for every cell
  // up front, because the neighbor lookup below needs the numbering of cells that come later in
  // the loop.
#ifndef NDEBUG
  const auto recomputedOrders = canonicalVertexOrders(meshTopology, meshGeometry);
#endif
  std::vector<std::array<std::uint8_t, Cell::NumFaces>> pumlFaceMaps(cells.size());
  for (std::size_t i = 0; i < cells.size(); i++) {
    VertexOrder vertexOrder{};
    std::copy_n(vertexOrders.entity(i), Cell::NumVertices, vertexOrder.begin());
    assert(vertexOrder == recomputedOrders[i]);
    pumlFaceMaps[i] = canonicalFaceMap(vertexOrder);
  }

  // Compute everything local
  elements_.resize(cells.size());
  for (std::size_t i = 0; i < cells.size(); i++) {
    elements_[i].globalId = cellIdsAsInFile[i];
    elements_[i].localId = i;
    elements_[i].clusterId = clusterIds[i];
    elements_[i].timestep = timestep[i];

    const auto& pumlToSeisSol = pumlFaceMaps[i];
    const auto vertexOrder = vertexOrderFromFaceMap(pumlToSeisSol);

    // Vertices
    std::array<PUML::LocalId, Cell::NumVertices> geomVertices{};
    PUML::Downward::vertices(meshGeometry, cellsGeometry[i], geomVertices.data());
    applyVertexOrder(geomVertices, vertexOrder);
    std::copy(geomVertices.begin(), geomVertices.end(), elements_[i].vertices.begin());

    std::array<PUML::LocalId, Cell::NumVertices> topoVertices{};
    PUML::Downward::vertices(meshTopology, cells[i], topoVertices.data());
    applyVertexOrder(topoVertices, vertexOrder);

    // Neighbor information
    std::array<PUML::LocalId, Cell::NumFaces> faceids{};
    PUML::Downward::faces(meshTopology, cells[i], faceids.data());
    std::array<PUML::LocalId, Cell::NumFaces> neighbors{};
    PUML::Neighbor::face(meshTopology, i, neighbors.data());

    for (std::size_t j = 0; j < Cell::NumFaces; j++) {
      const auto faceTag = boundaryTags.tag(i, j);
      const bool isLocallyCorrect = checkMeshCorrectnessLocally<PumlTopology>(
          faces[faceids[j]], neighbors, j, faceTag, cellIdsAsInFile[i], faceMap);
      isMeshCorrect &= isLocallyCorrect;
      const auto side = pumlToSeisSol[j];

      if (neighbors[j] == PUML::InvalidLocalId) {
        elements_[i].neighbors[side] = OptionalSize();

        if (!faces[faceids[j]].isShared()) {
          // Boundary sides
          elements_[i].neighborRanks[side] = rank;
        } else {
          // MPI Boundary
          neighborInfo[faces[faceids[j]].shared()[0]].push_back(faceids[j]);

          elements_[i].neighborRanks[side] = faces[faceids[j]].shared()[0];
        }
      } else {
        logassert(static_cast<std::size_t>(neighbors[j]) < cells.size());

        elements_[i].neighbors[side] = neighbors[j];

        std::array<PUML::LocalId, Cell::NumFaces> nfaces{};
        PUML::Neighbor::face(meshTopology, neighbors[j], nfaces.data());
        const auto* back = std::find(nfaces.begin(), nfaces.end(), i);
        logassert(back != nfaces.end());

        const auto& neighborPumlToSeisSol = pumlFaceMaps[neighbors[j]];
        elements_[i].neighborSides[side] = neighborPumlToSeisSol[back - nfaces.begin()];

        const auto firstVertex = topoVertices[FirstFaceVertex[side]];

        std::array<PUML::LocalId, Cell::NumVertices> nvertices{};
        PUML::Downward::vertices(meshTopology, cells[neighbors[j]], nvertices.data());
        applyVertexOrder(nvertices, vertexOrderFromFaceMap(neighborPumlToSeisSol));
        const auto* neighborFirstVertex =
            std::find(nvertices.begin(), nvertices.end(), firstVertex);
        logassert(neighborFirstVertex != nvertices.end());

        elements_[i].sideOrientations[side] =
            FaceVertexToOrientation[elements_[i].neighborSides[side]]
                                   [neighborFirstVertex - nvertices.begin()];
        logassert(elements_[i].sideOrientations[side] == 0);

        elements_[i].neighborRanks[side] = rank;
      }

      const auto bcCurrentFace = faceMap.at(faceTag);

      elements_[i].boundaries[side] = bcCurrentFace.value_or(FaceType::Regular);
      elements_[i].faultTags[side] = faceTag;
      elements_[i].mpiIndices[side] = 0;
    }

    elements_[i].group = group[i];
  }
  if (!isMeshCorrect) {
    logError() << "Found at least one broken face in the mesh, see errors above for a more "
                  "detailled analysis.";
  }

  // Exchange ghost layer information and generate neighbor list
  std::vector<std::vector<char>> copySide(neighborInfo.size());
  std::vector<std::vector<char>> ghostSide(neighborInfo.size());
  std::vector<std::vector<unsigned long>> copyFirstVertex(neighborInfo.size());
  std::vector<std::vector<unsigned long>> ghostFirstVertex(neighborInfo.size());

  std::vector<MPI_Request> requests(neighborInfo.size() * 4);

  std::unordered_set<PUML::LocalId> t;
  std::size_t sum = 0;

  for (auto [k, info] : seissol::common::enumerate(neighborInfo)) {
    // Need to sort the neighborInfo vectors once
    std::sort(info.second.begin(), info.second.end(), [&](auto a, auto b) {
      return meshTopology.faces()[a].gid() < meshTopology.faces()[b].gid();
    });

    t.insert(info.second.begin(), info.second.end());
    sum += info.second.size();

    // Create MPI neighbor list
    addMPINeighor(meshTopology, info.first, info.second, pumlFaceMaps);

    copySide[k].resize(info.second.size());
    ghostSide[k].resize(info.second.size());
    copyFirstVertex[k].resize(info.second.size());
    ghostFirstVertex[k].resize(info.second.size());

    MPI_Irecv(ghostSide[k].data(),
              info.second.size(),
              MPI_CHAR,
              info.first,
              0,
              Mpi::mpi.comm(),
              &requests[k]);
    MPI_Irecv(ghostFirstVertex[k].data(),
              info.second.size(),
              MPI_UNSIGNED_LONG,
              info.first,
              0,
              Mpi::mpi.comm(),
              &requests[neighborInfo.size() + k]);

    // Neighbor side
    for (std::size_t i = 0; i < info.second.size(); i++) {
      // The side of boundary
      std::array<PUML::LocalId, 2> cellIds{};
      PUML::Upward::cells(meshTopology, faces[info.second[i]], cellIds.data());
      const auto pumlSide =
          PUML::Downward::faceSide(meshTopology, cells[cellIds[0]], info.second[i]);
      logassert(pumlSide >= 0 && static_cast<std::size_t>(pumlSide) < Cell::NumFaces);

      // Send the SeisSol side under the canonical numbering: the receiving rank cannot translate
      // a PUML side, because the numbering is per cell now.
      const auto vertexOrder = vertexOrderFromFaceMap(pumlFaceMaps[cellIds[0]]);
      const auto side = pumlFaceMaps[cellIds[0]][pumlSide];
      copySide[k][i] = static_cast<char>(side);

      std::array<PUML::LocalId, Cell::NumVertices> topoVertices{};
      PUML::Downward::vertices(meshTopology, cells[cellIds[0]], topoVertices.data());
      applyVertexOrder(topoVertices, vertexOrder);

      // First vertex of the face on the boundary
      const auto firstVertex = topoVertices[FirstFaceVertex[side]];
      copyFirstVertex[k][i] = vertices[firstVertex].gid();

      // Set the MPI index
      logassert(elements_[cellIds[0]].mpiIndices[side] == 0);
      elements_[cellIds[0]].mpiIndices[side] = i;
    }

    MPI_Isend(copySide[k].data(),
              info.second.size(),
              MPI_CHAR,
              info.first,
              0,
              Mpi::mpi.comm(),
              &requests[neighborInfo.size() * 2 + k]);
    MPI_Isend(copyFirstVertex[k].data(),
              info.second.size(),
              MPI_UNSIGNED_LONG,
              info.first,
              0,
              Mpi::mpi.comm(),
              &requests[neighborInfo.size() * 3 + k]);
  }
  logassert(t.size() == sum);

  MPI_Waitall(neighborInfo.size() * 4, requests.data(), MPI_STATUSES_IGNORE);

  for (auto [k, info] : seissol::common::enumerate(neighborInfo)) {
    for (std::size_t i = 0; i < info.second.size(); i++) {
      // Set neighbor side
      std::array<PUML::LocalId, 2> cellIds{};
      PUML::Upward::cells(meshTopology, faces[info.second[i]], cellIds.data());
      logassert(cellIds[1] == PUML::InvalidLocalId);

      // the linters demanded a double cast here; both sides are already SeisSol sides under the
      // canonical numbering, so no translation is applied to either of them
      const auto side = static_cast<std::size_t>(static_cast<unsigned char>(copySide[k][i]));
      const auto gSide = static_cast<std::size_t>(static_cast<unsigned char>(ghostSide[k][i]));
      elements_[cellIds[0]].neighborSides[side] = gSide;

      // Set side sideOrientation
      std::array<unsigned long, Cell::NumVertices> nvertices{};
      PUML::Downward::gvertices(meshTopology, cells[cellIds[0]], nvertices.data());
      applyVertexOrder(nvertices, vertexOrderFromFaceMap(pumlFaceMaps[cellIds[0]]));

      const auto* localFirstVertex =
          std::find(nvertices.begin(), nvertices.end(), ghostFirstVertex[k][i]);
      logassert(localFirstVertex != nvertices.end());

      elements_[cellIds[0]].sideOrientations[side] =
          FaceVertexToOrientation[side][localFirstVertex - nvertices.begin()];
      logassert(elements_[cellIds[0]].sideOrientations[side] == 0);
    }
  }

  // Set vertices
  vertices_.resize(verticesGeometry.size());
  std::vector<PUML::LocalId> preElements;
  for (std::size_t i = 0; i < verticesGeometry.size(); i++) {
    std::copy_n(
        verticesGeometry[i].coordinate(), vertices_[i].coords.size(), vertices_[i].coords.begin());

    PUML::Upward::cells(meshGeometry, verticesGeometry[i], preElements);
    vertices_[i].elements.resize(preElements.size());
    std::copy(preElements.begin(), preElements.end(), vertices_[i].elements.begin());
  }

  // the neighborSide needs to be _inferred_ here.
  for (auto& [_, neighbor] : mpiNeighbors_) {
    for (auto& element : neighbor.elements) {
      element.neighborSide = elements_[element.localElement].neighborSides[element.localSide];
    }
  }
}

void PUMLReader::addMPINeighor(
    const PumlMesh& meshTopology,
    int rank,
    const std::vector<PUML::LocalId>& faces,
    const std::vector<std::array<std::uint8_t, Cell::NumFaces>>& pumlFaceMaps) {
  const std::size_t id = mpiNeighbors_.size();
  MPINeighbor& neighbor = mpiNeighbors_[rank];

  neighbor.localID = id;
  neighbor.elements.resize(faces.size());

  for (std::size_t i = 0; i < faces.size(); i++) {
    std::array<PUML::LocalId, 2> cellIds{};
    PUML::Upward::cells(meshTopology, meshTopology.faces()[faces[i]], cellIds.data());

    neighbor.elements[i].localElement = cellIds[0];

    neighbor.elements[i].neighborElement = i;

    std::array<PUML::LocalId, Cell::NumFaces> sides{};
    PUML::Downward::faces(meshTopology, meshTopology.cells()[cellIds[0]], sides.data());
    neighbor.elements[i].localSide = [&]() -> std::size_t {
      for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
        if (sides[f] == faces[i]) {
          return pumlFaceMaps[cellIds[0]][f];
        }
      }
      throw;
    }();
  }
}

bool PUMLReader::inlineTimestepCompute() const { return true; }

bool PUMLReader::inlineClusterCompute() const { return true; }

} // namespace seissol::geometry
