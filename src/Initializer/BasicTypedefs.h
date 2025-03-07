// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_BASICTYPEDEFS_H_
#define SEISSOL_SRC_INITIALIZER_BASICTYPEDEFS_H_

namespace seissol {

constexpr int DataTagOffset = 2;

enum class HaloType { Ghost, Copy, Interior };

enum class TimeClustering {
  // global time stepping
  Single = 0,
  // online clustering resulting in a multi-rate scheme
  MultiRate = 2,
};

// face types
// Note: When introducting new types also change
// int seissol::initializer::time_stepping::LtsWeights::getBoundaryCondition
// and PUMLReader. Otherwise it might become a DR face...
enum class FaceType {
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

  // periodic boundary
  Periodic = 6,

  // analytical boundary (from initial cond.)
  Analytical = 7
};

// Once the FaceType enum is updated, make sure to update these methods here as well.

// Checks if a face type is an internal face (i.e. there are two cells adjacent to it).
// That includes all interior and dynamic rupture faces, but also periodic faces.
constexpr bool isInternalFaceType(FaceType faceType) {
  return faceType == FaceType::Regular || faceType == FaceType::DynamicRupture ||
         faceType == FaceType::Periodic;
}

// Checks if a face type is an external boundary face (i.e. there is only one cell adjacent to it).
constexpr bool isExternalBoundaryFaceType(FaceType faceType) {
  return faceType == FaceType::FreeSurface || faceType == FaceType::FreeSurfaceGravity ||
         faceType == FaceType::Dirichlet || faceType == FaceType::Analytical;
}

enum class ComputeGraphType {
  AccumulatedVelocities = 0,
  StreamedVelocities,
  NeighborIntegral,
  DynamicRuptureInterface,
  Count
};

} // namespace seissol

#endif // SEISSOL_SRC_INITIALIZER_BASICTYPEDEFS_H_
