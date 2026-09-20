// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_BASICTYPEDEFS_H_
#define SEISSOL_SRC_INITIALIZER_BASICTYPEDEFS_H_

#include <array>
#include <cstdint>
#include <string_view>
namespace seissol {

enum class HaloType { Ghost, Copy, Interior };

/*
  Face types that a tetrahedron face may assume.

  Note: When introducing new types, also add them to FaceTypes and faceTypeName below, and
  extend the classification methods further down. Everything that enumerates or names face
  types derives from those. Additionally, every material model has to state whether the new
  type is defined for it (MaterialSetup<T>::supportsFaceType), and every solver whether it
  implements it (Solver::implementsFaceType).
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

// All face types. The enum values are not contiguous, so iterate this array instead of a
// numeric range.
constexpr std::array<FaceType, 7> FaceTypes = {FaceType::Regular,
                                               FaceType::FreeSurface,
                                               FaceType::FreeSurfaceGravity,
                                               FaceType::DynamicRupture,
                                               FaceType::Dirichlet,
                                               FaceType::Outflow,
                                               FaceType::Analytical};

// The name a face type is referred to by in the mesh face map and in diagnostics.
constexpr std::string_view faceTypeName(FaceType faceType) {
  switch (faceType) {
  case FaceType::Regular:
    return "regular";
  case FaceType::FreeSurface:
    return "freeSurface";
  case FaceType::FreeSurfaceGravity:
    return "freeSurfaceGravity";
  case FaceType::DynamicRupture:
    return "dynamicRupture";
  case FaceType::Dirichlet:
    return "dirichlet";
  case FaceType::Outflow:
    return "outflow";
  case FaceType::Analytical:
    return "analytical";
  }
  return "unknown";
}

/**
 * Whether a face type is available in a given context, and if not, why.
 *
 * Two independent contexts use this: a material model states which boundary conditions are
 * defined for it, and a solver states which ones its kernels implement. The distinction
 * matters because the two are fixed by different things -- a missing solver implementation
 * can be added, while a boundary condition that is not formulated for a material cannot.
 */
struct FaceTypeSupport {
  bool supported{true};
  std::string_view reason;
};

constexpr FaceTypeSupport faceTypeSupported() { return {}; }

constexpr FaceTypeSupport faceTypeUnsupported(std::string_view reason) { return {false, reason}; }

enum class BCType {
  // an internal face, with a neighbor
  Internal,

  // a boundary face with a fake neighbor (it needs the Neighbor kernel)
  ExternalFake,

  // a boundary face without a fake neighbor (can handle everything in the Local kernel)
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
      faceType == FaceType::Dirichlet || faceType == FaceType::Analytical ||
      faceType == FaceType::Outflow) {
    return BCType::ExternalNone;
  }

  // currently, there is no ExternalFake BC

  // should never happen, unless you forgot to implement something
  return BCType::Unknown;
}

// Checks if a face type is an internal face (i.e. there are two cells adjacent to it).
// That includes all interior and dynamic rupture faces, but also periodic faces.
constexpr bool isInternalFaceType(FaceType faceType) {
  return getBCType(faceType) == BCType::Internal;
}

// Checks if a face type builds its boundary state in the nodal face basis and applies it
// through the local kernel, as opposed to having it folded into the flux solver matrices
// during setup.
constexpr bool requiresNodalFlux(FaceType faceType) {
  return faceType == FaceType::FreeSurfaceGravity || faceType == FaceType::Dirichlet ||
         faceType == FaceType::Analytical;
}

// Checks if a face type has to be evaluated with the Godunov flux, regardless of the
// numerical flux selected in the parameter file.
constexpr bool enforcesGodunovFlux(FaceType faceType) {
  return faceType == FaceType::FreeSurface || faceType == FaceType::FreeSurfaceGravity ||
         faceType == FaceType::Analytical || faceType == FaceType::Outflow;
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
