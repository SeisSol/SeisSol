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

using GlobalElemId = size_t;
using LocalElemId = int;   // TODO(David): size_t, maybe, once the Netcdf Reader is gone
using LocalVertexId = int; // TODO(David): size_t, maybe, once the Netcdf Reader is gone
using SideId = int;        // TODO(David): int8_t , once the Netcdf Reader is gone

using ElemVertices = LocalVertexId[Cell::NumVertices];
using ElemNeighbors = LocalElemId[Cell::NumFaces];
using ElemNeighborSides = SideId[Cell::NumFaces];
using ElemSideOrientations = int[Cell::NumFaces];
using ElemBoundaries = int[Cell::NumFaces];
using ElemNeighborRanks = int[Cell::NumFaces]; // type prescribed by MPI
/** The index of this element (side) in the communication array */
using ElemMPIIndices = int[Cell::NumFaces];
using ElemGroup = int;
using ElemFaultTags = int[Cell::NumFaces];

using VrtxCoords = double[Cell::Dim];

struct Element {
  GlobalElemId globalId;
  LocalElemId localId;
  ElemVertices vertices;
  ElemNeighbors neighbors;
  ElemNeighborSides neighborSides;
  ElemSideOrientations sideOrientations;
  /** Domain boundary condition, or 0 for inner elements and MPI boundaries */
  ElemBoundaries boundaries;
  ElemNeighborRanks neighborRanks;
  ElemMPIIndices mpiIndices;
  ElemMPIIndices mpiFaultIndices;
  /** Material of the element */
  ElemGroup group;
  ElemFaultTags faultTags; // member of struct Element
  int clusterId;
  double timestep;
};

struct Vertex {
  VrtxCoords coords{};
  /** Elements sharing this neighbor */
  std::vector<LocalElemId> elements;
};

struct MPINeighborElement {
  /** Local number of the local element */
  LocalElemId localElement;
  /** Side of the local element */
  SideId localSide;
  /** Global number neighbor element */
  LocalElemId neighborElement;
  /** Side of the neighbor element */
  SideId neighborSide;
};

struct Fault {
  GlobalElemId globalId;
  /** The element which contains this fault */
  LocalElemId element;
  /** The side of the element */
  SideId side;

  GlobalElemId neighborGlobalId;
  LocalElemId neighborElement;
  SideId neighborSide;
  int tag;

  /** Normal of the fault face */
  VrtxCoords normal;

  /** Rotation matrix */
  VrtxCoords tangent1;
  VrtxCoords tangent2;
};

struct MPINeighbor {
  /** Local ID of the MPI neighbor */
  LocalElemId localID{};

  std::vector<MPINeighborElement> elements;
};

} // namespace seissol

#endif // SEISSOL_SRC_GEOMETRY_MESHDEFINITION_H_
