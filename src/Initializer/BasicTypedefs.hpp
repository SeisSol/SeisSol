#ifndef SEISSOL_BASICTYPEDEFS_HPP
#define SEISSOL_BASICTYPEDEFS_HPP

#include <Kernels/precision.hpp>

enum mpiTag { localIntegrationData = 0, neighboringIntegrationData = 1, timeData = 2 };

enum TimeClustering {
  // global time stepping
  single = 0,
  // offline clustering computed in pre-processing
  offline = 1,
  // online clustering resulting in a multi-rate scheme
  multiRate = 2,
  // online clustering aiming at LTS for slithers only
  slithers = 3
};

// face types
// Note: When introducting new types also change
// int seissol::initializers::time_stepping::LtsWeights::getBoundaryCondition
// and PUMLReader. Otherwise it might become a DR face...
enum class FaceType {
  // regular: inside the computational domain
  regular = 0,

  // free surface boundary
  freeSurface = 1,

  // free surface boundary with gravity
  freeSurfaceGravity = 2,

  // dynamic rupture boundary
  dynamicRupture = 3,

  // Dirichlet boundary
  dirichlet = 4,

  // absorbing/outflow boundary
  outflow = 5,

  // periodic boundary
  periodic = 6,

  // analytical boundary (from initial cond.)
  analytical = 7
};

// Once the FaceType enum is updated, make sure to update these methods here as well.

// Checks if a face type is an internal face (i.e. there are two cells adjacent to it).
// That includes all interior and dynamic rupture faces, but also periodic faces.
constexpr bool isInternalFaceType(FaceType faceType) {
  return faceType == FaceType::regular || faceType == FaceType::dynamicRupture ||
         faceType == FaceType::periodic;
}

// Checks if a face type is an external boundary face (i.e. there is only one cell adjacent to it).
constexpr bool isExternalBoundaryFaceType(FaceType faceType) {
  return faceType == FaceType::freeSurface || faceType == FaceType::freeSurfaceGravity ||
         faceType == FaceType::dirichlet || faceType == FaceType::analytical;
}

enum SystemType { Host = 0, Device = 1 };

enum class ComputeGraphType {
  LocalIntegral = 0,
  AccumulatedVelocities,
  StreamedVelocities,
  NeighborIntegral,
  DynamicRuptureInterface,
  Count
};

#endif // SEISSOL_BASICTYPEDEFS_HPP
