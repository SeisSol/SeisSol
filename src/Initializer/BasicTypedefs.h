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

// face types
// Note: When introducting new types also change
// int seissol::initializer::time_stepping::LtsWeights::getBoundaryCondition
// and PUMLReader. Otherwise it might become a DR face...
enum class FaceType : uint8_t {
  // regular: inside the computational domain
  Regular = 0,

  // free surface boundary
  FreeSurface = 1,

  // free surface boundary with gravity
  FreeSurfaceGravity = 2,

  // dynamic rupture boundary
  DynamicRupture = 3,

  // Dirichlet boundary
  Dirichlet = 4,

  // absorbing/outflow boundary
  Outflow = 5,

  // periodic boundary (obsolete; now equivalent to Regular)
  // Periodic = 6,

  // analytical boundary (from initial cond.)
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
