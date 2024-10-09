/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de,
 *http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013-2014, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Defines mesh data structures
 **/

#ifndef MESH_DEFINITION_H
#define MESH_DEFINITION_H

#include <cstddef>
#include <vector>

namespace seissol {

using GlobalElemId = size_t;
using LocalElemId = int;   // TODO(David): size_t, maybe, once the Netcdf Reader is gone
using LocalVertexId = int; // TODO(David): size_t, maybe, once the Netcdf Reader is gone
using SideId = int;        // TODO(David): int8_t , once the Netcdf Reader is gone

using ElemVertices = LocalVertexId[4];
using ElemNeighbors = LocalElemId[4];
using ElemNeighborSides = SideId[4];
using ElemSideOrientations = int[4];
using ElemBoundaries = int[4];
using ElemNeighborRanks = int[4]; // type prescribed by MPI
/** The index of this element (side) in the communication array */
using ElemMPIIndices = int[4];
using ElemGroup = int;
using ElemFaultTags = int[4];

using VrtxCoords = double[3];

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
};

struct Vertex {
  VrtxCoords coords;
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
  LocalElemId localID;

  std::vector<MPINeighborElement> elements;
};

} // namespace seissol

#endif // MESH_DEFINITION_H
