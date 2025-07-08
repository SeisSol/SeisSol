// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_GEOMETRY_MESHDEFINITION_H_
#define SEISSOL_SRC_GEOMETRY_MESHDEFINITION_H_

#include <Common/Constants.h>
#include <cstddef>
#include <vector>

namespace seissol {

using CoordinateT = std::array<double, Cell::Dim>;

struct Element {
  // global cell ID
  size_t globalId{};
  size_t localId{};
  std::array<size_t, Cell::NumVertices> vertices{};
  std::array<size_t, Cell::NumFaces> neighbors{};
  std::array<int8_t, Cell::NumFaces> neighborSides{};
  std::array<int8_t, Cell::NumFaces> sideOrientations{};
  /** Domain boundary condition, or 0 for inner elements and MPI boundaries */
  std::array<int, Cell::NumFaces> boundaries{};
  std::array<int, Cell::NumFaces> neighborRanks{};
  std::array<size_t, Cell::NumFaces> mpiIndices{};
  std::array<size_t, Cell::NumFaces> mpiFaultIndices{};
  /** Material of the element */
  int group{};
  std::array<int, Cell::NumFaces> faultTags{};
  int clusterId{};
  double timestep{};
};

struct Vertex {
  CoordinateT coords{};
  /** Elements sharing this neighbor */
  std::vector<size_t> elements;
};

struct MPINeighborElement {
  /** Local number of the local element */
  size_t localElement{};
  /** Side of the local element */
  int8_t localSide{};
  /** Global number neighbor element */
  size_t neighborElement{};
  /** Side of the neighbor element */
  int8_t neighborSide{};
};

struct Fault {
  size_t globalId{};
  /** The element which contains this fault */
  size_t element{};
  /** The side of the element */
  int8_t side{};

  size_t neighborGlobalId{};
  size_t neighborElement{};
  int8_t neighborSide{};
  int tag{};

  /** Normal of the fault face */
  CoordinateT normal{};

  /** Rotation matrix */
  CoordinateT tangent1{};
  CoordinateT tangent2{};
};

struct MPINeighbor {
  /** Local ID of the MPI neighbor */
  size_t localID{};

  std::vector<MPINeighborElement> elements;
};

} // namespace seissol

#endif // SEISSOL_SRC_GEOMETRY_MESHDEFINITION_H_
