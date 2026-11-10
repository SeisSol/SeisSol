// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PREFLIGHTMESH_H_
#define SEISSOL_SRC_INITIALIZER_PREFLIGHTMESH_H_

#include "Geometry/MeshDefinition.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/FaceMap.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Initializer/Parameters/MeshParameters.h"

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace seissol::initializer::meshscan {

/**
 * @brief Points that the run locates in the cells of the mesh, as findUniqueMeshIds does: the
 * off-fault receivers or the point sources.
 */
struct CellPoints {
  std::vector<Eigen::Vector3d> points;
  //! the tolerance the run locates them with
  double tolerance{0};

  //! What the scan finds out about a point.
  struct Found {
    //! the cell the run takes for the point, by its index in the mesh file, if there is one
    std::optional<std::size_t> cell;
    //! the height of the free surface above or below the point where there is any, i.e. the
    //! highest point of the faces of the free surface that the vertical through the point meets
    std::optional<double> surface;
  };
  //! filled in by scanMesh, one entry per point
  std::vector<Found> found;
};

/**
 * @brief The on-fault receivers, which the run puts onto the fault face whose prism along its
 * normal holds them, the closest one if there are several (see PickPointBuilder).
 */
struct FaultPoints {
  std::vector<CoordinateT> points;

  //! Where the run puts a point.
  struct Found {
    //! the cell of the plus side of the fault face, by its index in the mesh file
    std::size_t cell{0};
    //! the distance from the point to the face, by which the run moves it onto the face
    double distance{0};
    //! the square root of the area of the face, i.e. about the length of its edges
    double faceSize{0};
  };
  //! filled in by scanMesh, one entry per point; empty where the run finds no face
  std::vector<std::optional<Found>> found;
};

//! @brief The cells of a chunk of the mesh file, as the run takes them.
struct CellChunk {
  //! the index of the first of them in the file
  std::size_t first{0};
  //! their vertices, displaced and scaled, in the order the run puts them in (zero for a cell
  //! with vertices the file does not have)
  std::vector<std::array<std::array<double, 3>, 4>> vertices;
  std::vector<int> groups;
  //! whether the file has the vertices of the cell
  std::vector<std::uint8_t> valid;
};

//! @brief What the scan needs to know of the run.
struct MeshScanSettings {
  std::string file;
  //! Auto takes the one the file gives (see boundaryFormatOf)
  parameters::BoundaryFormat boundaryFormat{parameters::BoundaryFormat::Auto};
  //! the mesh is displaced first and scaled then, as the run does (see MeshParameters)
  Eigen::Vector3d displacement{Eigen::Vector3d::Zero()};
  Eigen::Matrix3d scaling{Eigen::Matrix3d::Identity()};
  FaceMap faceMap;
  //! whether the faces of dynamic rupture are such; the run takes them for regular ones else
  bool dynamicRupture{false};
  //! how the run tells the sides of a fault face (see MeshReader::extractFaultInformation)
  Eigen::Vector3d referencePoint{Eigen::Vector3d::Zero()};
  parameters::RefPointMethod refPointMethod{parameters::RefPointMethod::Point};
  //! the number of cells read at once
  std::size_t chunk{std::size_t{1} << 20};
  //! called with every chunk of cells, one after the other, if given
  std::function<void(const CellChunk&)> visitCells;
};

//! @brief What the scan finds out about the mesh as a whole.
struct MeshScan {
  std::size_t cells{0};
  std::size_t vertices{0};
  parameters::BoundaryFormat boundaryFormat{parameters::BoundaryFormat::Auto};
  //! the bounding box of the vertices, displaced and scaled
  Eigen::Vector3d min{Eigen::Vector3d::Zero()};
  Eigen::Vector3d max{Eigen::Vector3d::Zero()};

  //! the faces by their boundary condition tag, where it is not zero; a face between two cells
  //! counts for both of them
  std::map<std::uint32_t, std::size_t> tags;
  //! the tags that the face map knows no face type for, and how many faces have them
  std::map<std::uint32_t, std::size_t> unknownTags;
  //! the faces by their face type, as the run takes them, but the regular ones; counted as the tags
  std::map<FaceType, std::size_t> faceTypes;
  //! the cells by their group
  std::map<int, std::size_t> groups;
  //! the face types the cells of each group have, as masks of internal::faceTypeBit
  std::map<int, std::uint32_t> groupFaceTypes;
  //! the faces of dynamic rupture, each counted once
  std::size_t faultFaces{0};

  //! the cells whose vertices index beyond the vertices of the file
  std::size_t invalidCells{0};
  //! the cells the run rejects for their orientation (see MeshReader::verifyMeshOrientation),
  //! e.g. as they have no volume; the first of them, by its index in the file
  std::size_t misorientedCells{0};
  std::optional<std::size_t> firstMisorientedCell;
  //! the radius of the largest sphere in the smallest cell, and that cell
  double smallestInsphere{0};
  std::size_t smallestCell{0};
  double largestInsphere{0};

  //! the memory the scan held, mostly the vertices
  std::size_t bytes{0};
  double seconds{0};
};

/**
 * @brief Reads the PUML mesh file once, cell by cell, and takes what the run would find out about
 * it only after reading and partitioning it: the boundary conditions and groups of the cells, cells
 * that the run rejects, and where the points lie.
 *
 * Runs on one process, without MPI, with the threads of OpenMP. It holds the vertices of the mesh
 * and a chunk of the cells only, not the mesh.
 *
 * The points are located by the same tests as in the run (CellPlanes, isInsideFace), on the cells
 * in the vertex order the run puts them in (canonicalVertexOrder).
 *
 * Logs an error, i.e. stops, if the file cannot be read.
 */
MeshScan scanMesh(const MeshScanSettings& settings,
                  CellPoints& receivers,
                  CellPoints& sources,
                  FaultPoints& faultReceivers);

} // namespace seissol::initializer::meshscan

#endif // SEISSOL_SRC_INITIALIZER_PREFLIGHTMESH_H_
