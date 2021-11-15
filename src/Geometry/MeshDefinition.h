/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
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

#include <vector>

typedef int ElemVertices[4];
typedef int ElemNeighbors[4];
typedef int ElemNeighborSides[4];
typedef int ElemSideOrientations[4];
typedef int ElemBoundaries[4];
typedef int ElemNeighborRanks[4];
/** The index of this element (side) in the communication array */
typedef int ElemMPIIndices[4];
typedef int ElemMaterial;

typedef int ElemFaultTags[4];

struct Element {
	int localId;
	ElemVertices vertices;
	int rank;
	ElemNeighbors neighbors;
	ElemNeighborSides neighborSides;
	ElemSideOrientations sideOrientations;
	/** Domain boundary condition, or 0 for inner elements and MPI boundaries */
	ElemBoundaries boundaries;
	ElemNeighborRanks neighborRanks;
	ElemMPIIndices mpiIndices;
	ElemMPIIndices mpiFaultIndices;
	/** Material of the element */
	ElemMaterial material;
   ElemFaultTags faultTags; // member of struct Element
};

typedef double VrtxCoords[3];

struct Vertex {
	VrtxCoords coords;
	/** Elements sharing this neighbor */
	std::vector<int> elements;
};

struct MPINeighborElement {
	/** Local number of the local element */
	int localElement;
	/** Side of the local element */
	int localSide;
	/** Global number neighbor element */
	int neighborElement;
	/** Side of the neighbor element */
	int neighborSide;

	/**
	 * Sort elements by according to the local order
	 *
	 * @todo Remove this function, since we need to sort by local and neighbor ids
	 */
	bool operator<(const MPINeighborElement &other) const
	{
		return (localElement < other.localElement)
				|| (localElement == other.localElement
						&& localSide < other.localSide);
	}
};

struct Fault {
	/** The element which contains this fault */
	int element;
	/** The side of the element */
	int side;

	int neighborElement;
	int neighborSide;

	/** Normal of the fault face */
	VrtxCoords normal;

	/** Rotation matrix */
	VrtxCoords tangent1;
	VrtxCoords tangent2;
};

struct MPINeighbor {
	/** Local ID of the MPI neighbor */
	int localID;

	std::vector<MPINeighborElement> elements;
};

#endif // MESH_DEFINITION_H
