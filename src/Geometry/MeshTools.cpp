/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2014, SeisSol Group
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
 **/

#include "MeshTools.h"

#include <algorithm>
#include <iostream>
#include <iterator>

#include <Eigen/Dense>

const int MeshTools::FACE2NODES[4][3] = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};
const int MeshTools::FACE2MISSINGNODE[4] = {3, 2, 1, 0};
const int MeshTools::NEIGHBORFACENODE2LOCAL[3] = {0, 2, 1};

void Plane::transform(const VrtxCoords globalCoords,
		      VrtxCoords localCoords) const {
  Eigen::Vector3d coordsVec;
  coordsVec << globalCoords[0], globalCoords[1], globalCoords[2];
  Eigen::Matrix3d invMapping;
  invMapping << normal[0], tangent1[0], tangent2[0],
    normal[1], tangent1[1], tangent2[1],
    normal[2], tangent1[2], tangent2[2];
  Eigen::Vector3d x = invMapping.colPivHouseholderQr().solve(coordsVec);
  localCoords[0] = x[0];
  localCoords[1] = x[1];
  localCoords[2] = x[2];
}

void Plane::transformBack(const VrtxCoords localCoords,
			  VrtxCoords globalCoords) const {
  Eigen::Matrix3d invMapping;
  invMapping << normal[0], tangent1[0], tangent2[0],
    normal[1], tangent1[1], tangent2[1],
    normal[2], tangent1[2], tangent2[2];
  Eigen::Vector3d coordsVec;
  coordsVec << localCoords[0], localCoords[1], localCoords[2];
  Eigen::Vector3d x = invMapping * coordsVec;
  globalCoords[0] = x[0];
  globalCoords[1] = x[1];
  globalCoords[2] = x[2];
}

PlaneFragment::PlaneFragment(const std::unordered_set<PeriodicVertex> &points,
			     Plane plane)
  : plane(plane) {
  // Find constant dimension and non constant dimensions.
  constDim = plane.getConstantDim();
  std::cout << plane.normal[constDim];
  int otherDims[2];
  plane.getOtherDims(otherDims);

  // Convert to 2D points
  auto pointsPlane = std::vector<PointPlane>{};
  for (const auto& point : points) {
    auto curPointPlane = PointPlane{};
    curPointPlane.coords[0] = point.globalCoords[otherDims[0]];
    curPointPlane.coords[1] = point.globalCoords[otherDims[1]];
    pointsPlane.push_back(curPointPlane);
  }
  convexHull = computeConvexHull2D(pointsPlane);

}

PlaneFragment::PlaneFragment() {

}

PlaneFragmentsSerialized::PlaneFragmentsSerialized(const std::vector<double> &dataBuffer,
						   const std::vector<int> &lengthsBuffer)
  : data(dataBuffer.data()),
    lengths(lengthsBuffer.data()),
    dataSize(dataBuffer.size()),
    lengthsSize(lengthsBuffer.size()){
}


void encodePlaneFragment(const PlaneFragment& planeFragment,
		  std::vector<double> &dataBuffer,
		  std::vector<int> &lengthsBuffer) {
  const auto &plane = planeFragment.plane;
  const auto &hull = planeFragment.convexHull;
  // Normal + Tangent1 + Tangent2 + Point + 2 * convexHull Points
  const auto fragSize = 3 + 3 + 3 + 3 + 2 * hull.size();
  dataBuffer.reserve(dataBuffer.size() + fragSize);

  auto bi = std::back_inserter(dataBuffer);
  std::copy_n(plane.normal, 3, bi);
  std::copy_n(plane.point, 3, bi);
  std::copy_n(plane.tangent1, 3, bi);
  std::copy_n(plane.tangent2, 3, bi);
  for (const auto &h : hull) {
    std::copy_n(h.coords, 2, bi);
  }

  lengthsBuffer.push_back(fragSize);

}

std::vector<PlaneFragment> decodePlaneFragments(const PlaneFragmentsSerialized &serialized,
						const std::vector<int> &planesPerNode) {
  auto planeFragments = std::vector<PlaneFragment>{};
  
  const auto &data = serialized.data;
  const auto &lengths = serialized.lengths;

  const auto dataSize = serialized.dataSize;

  // First decode the first one.
  auto lengthsIdx = 0;
  auto offset = 0;
  
  auto curRank = 0;
  auto curLocalPlaneIdx = 0;
  while (offset < serialized.dataSize) {
    assert(lengthsIdx < serialized.lengthsSize);
    PlaneFragment curFragment;
    // Decode plane.
    Plane plane;
    std::copy_n(&data[offset], 3, plane.normal);
    offset += 3;
    std::copy_n(&data[offset], 3, plane.point);
    offset += 3;
    std::copy_n(&data[offset], 3, plane.tangent1);
    offset += 3;
    std::copy_n(&data[offset], 3, plane.tangent2);
    offset += 3;
    curFragment.plane = std::move(plane);

    // Parse convex hull.
    const auto curHullSize = (lengths[lengthsIdx] - 4 * 3) / 2;
    curFragment.convexHull.resize(curHullSize);
    for (int i = 0; i < curHullSize; ++i) {
      std::copy_n(&data[offset], 2, curFragment.convexHull[i].coords);
      offset += 2;
    }
    curFragment.constDim = plane.getConstantDim();

    std::cout << "decode, offset, size, dataSize: " << offset << ", "
	      << serialized.lengths[lengthsIdx] << ", "
	      << serialized.dataSize << std::endl;

    ++lengthsIdx;

    if (planesPerNode[curRank] <= curLocalPlaneIdx) {
      curLocalPlaneIdx = 0;
      ++curRank;
    }
    std::cout << "Rank = " << curRank << " idx = " <<curLocalPlaneIdx
	      << " nplanes(rank) = " << planesPerNode[curRank]
	      << std::endl;
    curFragment.rank = curRank;
    curFragment.planeIdx = curLocalPlaneIdx;
    planeFragments.push_back(curFragment);
    ++curLocalPlaneIdx;

    assert(offset <= serialized.dataSize);
  }
  std::cout << "Found " << planeFragments.size() << " plane fragments" << std::endl;
    
  return planeFragments;
}

MPI_Datatype registerFaceMPI() {
  // MPI Types
  // Vertex
  int faceMpi_nitems = 4;
  int faceMpi_blocklengths[] = {3, 3, 3, 1};
  MPI_Datatype faceMpi_types[] = {MPI_DOUBLE,
				  MPI_DOUBLE,
				  MPI_DOUBLE,
				  MPI_INT};
  MPI_Datatype faceMpi_t;
  MPI_Aint faceMpi_offsets[] = {offsetof(FaceMPI, coords[0]),
				offsetof(FaceMPI, coords[1]),
				offsetof(FaceMPI, coords[2]),
				offsetof(FaceMPI, planeId)};

  MPI_Type_create_struct(faceMpi_nitems,
			 faceMpi_blocklengths,
			 faceMpi_offsets,
			 faceMpi_types,
			 &faceMpi_t);
  MPI_Type_commit(&faceMpi_t);
  return faceMpi_t;
}
