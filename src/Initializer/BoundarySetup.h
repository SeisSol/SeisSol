// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_BOUNDARYSETUP_H_
#define SEISSOL_SRC_INITIALIZER_BOUNDARYSETUP_H_

#include "Initializer/BasicTypedefs.h"

#include <string_view>

namespace seissol {

/**
 * What a face type needs from the setup, stated per face type in the way MaterialSetup states
 * what a material model needs. The answers are read while the cell-local matrices, the boundary
 * mappings and the storage are built; the kernels dispatch on the face type themselves.
 */
template <FaceType Type>
struct BoundarySetup {
  static constexpr std::string_view name() { return faceTypeName(Type); }

  static constexpr BCType bcType() { return getBCType(Type); }

  /// Whether the numerical flux is the Godunov one regardless of what the parameter file asks for.
  static constexpr bool enforcesGodunovFlux() { return false; }

  /// Whether the face carries its own setup data: the nodes, the rotation to the face-aligned
  /// basis, and whatever terms the condition is made of.
  static constexpr bool requiresFaceData() { return false; }

  /// Whether the ghost state the neighbouring flux solver acts on is already given in the
  /// face-aligned basis, so that the solver must not rotate it a second time.
  static constexpr bool usesFaceAlignedGhostState() { return false; }

  /// Whether the boundary state is built at the face nodes during the timestep and lifted back
  /// by the local kernel, as opposed to being known once the setup is done.
  static constexpr bool buildsNodalState() { return false; }

  /// Whether the condition is an affine map that is constant over the face, and therefore part
  /// of the local flux solver once the setup has folded it in.
  static constexpr bool foldsConditionIntoFluxSolver() { return false; }
};

template <>
struct BoundarySetup<FaceType::FreeSurface> : BoundarySetup<FaceType::Regular> {
  static constexpr std::string_view name() { return faceTypeName(FaceType::FreeSurface); }
  static constexpr BCType bcType() { return getBCType(FaceType::FreeSurface); }
  // the vanishing traction is built into the Godunov state itself
  static constexpr bool enforcesGodunovFlux() { return true; }
};

template <>
struct BoundarySetup<FaceType::FreeSurfaceGravity> : BoundarySetup<FaceType::Regular> {
  static constexpr std::string_view name() { return faceTypeName(FaceType::FreeSurfaceGravity); }
  static constexpr BCType bcType() { return getBCType(FaceType::FreeSurfaceGravity); }
  static constexpr bool enforcesGodunovFlux() { return true; }
  static constexpr bool requiresFaceData() { return true; }
  static constexpr bool usesFaceAlignedGhostState() { return true; }
  static constexpr bool foldsConditionIntoFluxSolver() { return true; }
};

template <>
struct BoundarySetup<FaceType::Dirichlet> : BoundarySetup<FaceType::Regular> {
  static constexpr std::string_view name() { return faceTypeName(FaceType::Dirichlet); }
  static constexpr BCType bcType() { return getBCType(FaceType::Dirichlet); }
  // the exterior state is prescribed, so the flux has to take the incoming characteristics
  // from it rather than blend it with the interior one
  static constexpr bool enforcesGodunovFlux() { return true; }
  static constexpr bool requiresFaceData() { return true; }
  static constexpr bool usesFaceAlignedGhostState() { return true; }
  static constexpr bool foldsConditionIntoFluxSolver() { return true; }
};

template <>
struct BoundarySetup<FaceType::Outflow> : BoundarySetup<FaceType::Regular> {
  static constexpr std::string_view name() { return faceTypeName(FaceType::Outflow); }
  static constexpr BCType bcType() { return getBCType(FaceType::Outflow); }
  // letting only the outgoing characteristics leave is what the Godunov state does
  static constexpr bool enforcesGodunovFlux() { return true; }
};

template <>
struct BoundarySetup<FaceType::Analytical> : BoundarySetup<FaceType::Regular> {
  static constexpr std::string_view name() { return faceTypeName(FaceType::Analytical); }
  static constexpr BCType bcType() { return getBCType(FaceType::Analytical); }
  static constexpr bool enforcesGodunovFlux() { return true; }
  static constexpr bool requiresFaceData() { return true; }
  static constexpr bool buildsNodalState() { return true; }
};

/// The answers of BoundarySetup for a face type that is only known at run time.
struct BoundaryProperties {
  std::string_view name;
  BCType bcType{BCType::Unknown};
  bool enforcesGodunovFlux{false};
  bool requiresFaceData{false};
  bool usesFaceAlignedGhostState{false};
  bool buildsNodalState{false};
  bool foldsConditionIntoFluxSolver{false};
};

namespace boundary_setup_detail {
template <FaceType Type>
constexpr BoundaryProperties properties() {
  return BoundaryProperties{BoundarySetup<Type>::name(),
                            BoundarySetup<Type>::bcType(),
                            BoundarySetup<Type>::enforcesGodunovFlux(),
                            BoundarySetup<Type>::requiresFaceData(),
                            BoundarySetup<Type>::usesFaceAlignedGhostState(),
                            BoundarySetup<Type>::buildsNodalState(),
                            BoundarySetup<Type>::foldsConditionIntoFluxSolver()};
}
} // namespace boundary_setup_detail

constexpr BoundaryProperties boundaryProperties(FaceType faceType) {
  switch (faceType) {
  case FaceType::Regular:
    return boundary_setup_detail::properties<FaceType::Regular>();
  case FaceType::FreeSurface:
    return boundary_setup_detail::properties<FaceType::FreeSurface>();
  case FaceType::FreeSurfaceGravity:
    return boundary_setup_detail::properties<FaceType::FreeSurfaceGravity>();
  case FaceType::DynamicRupture:
    return boundary_setup_detail::properties<FaceType::DynamicRupture>();
  case FaceType::Dirichlet:
    return boundary_setup_detail::properties<FaceType::Dirichlet>();
  case FaceType::Outflow:
    return boundary_setup_detail::properties<FaceType::Outflow>();
  case FaceType::Analytical:
    return boundary_setup_detail::properties<FaceType::Analytical>();
  }
  return BoundaryProperties{};
}

} // namespace seissol
#endif // SEISSOL_SRC_INITIALIZER_BOUNDARYSETUP_H_
