// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_BASICTYPEDEFS_H_
#define SEISSOL_SRC_INITIALIZER_BASICTYPEDEFS_H_

#include <cstdint>
namespace seissol {

enum class HaloType { Ghost, Copy, Interior };

/*
  Face types that a tetrahedron face may assume.

  Note: When introducting new types also change
  the defaultFaceMap() (FaceMap.cpp), the docs,
  and the methods below. Otherwise, it'll probably cause bugs and errors.
*/
enum class FaceType : uint8_t {
  // regular: inside the computational domain (interior face, linear)
  Regular = 0,

  // free surface boundary (boundary, linear)
  FreeSurface = 1,

  // free surface boundary with gravity (boundary, nonlinear)
  FreeSurfaceGravity = 2,

  // dynamic rupture boundary (interior face, nonlinear)
  DynamicRupture = 3,

  // Dirichlet boundary (boundary, nonlinear)
  Dirichlet = 4,

  // absorbing/outflow boundary (boundary, linear)
  Outflow = 5,

  // periodic boundary (obsolete; now equivalent to Regular)
  // Periodic = 6,

  // analytical boundary, taken from the initial conditions (boundary, nonlinear)
  Analytical = 7
};

enum class BCType {
  // an internal face, with a neighbor
  Internal,

  // a boundary face with a fake neighbor
  ExternalFake,

  // a boundary face without a fake neighbor
  ExternalNone,

  // unhandled face type
  Unknown
};

// Once the FaceType enum is updated, make sure to update these methods here as well.

constexpr BCType getBCType(FaceType faceType) {
  if (faceType == FaceType::Regular || faceType == FaceType::DynamicRupture) {
    return BCType::Internal;
  }
  if (faceType == FaceType::FreeSurface || faceType == FaceType::FreeSurfaceGravity ||
      faceType == FaceType::Dirichlet || faceType == FaceType::Analytical) {
    return BCType::ExternalFake;
  }
  if (faceType == FaceType::Outflow) {
    return BCType::ExternalNone;
  }

  // should never happen, unless you forgot to implement something
  return BCType::Unknown;
}

// Checks if a face type is an internal face (i.e. there are two cells adjacent to it).
// That includes all interior and dynamic rupture faces, but also periodic faces.
constexpr bool isInternalFaceType(FaceType faceType) {
  return getBCType(faceType) == BCType::Internal;
}

// Checks if a face type is an external boundary face (i.e. there is only one cell adjacent to it).
// (note that Outflow is purposefully excluded here)
constexpr bool isExternalBoundaryFaceType(FaceType faceType) {
  return getBCType(faceType) == BCType::ExternalFake;
}

enum class ComputeGraphType {
  AccumulatedVelocities = 0,
  StreamedVelocities,
  NeighborIntegral,
  DynamicRuptureInterface,
  Plasticity,
  Count
};

} // namespace seissol

#endif // SEISSOL_SRC_INITIALIZER_BASICTYPEDEFS_H_
