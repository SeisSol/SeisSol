// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_WAVEPROPAGATION_SOLVER_MATERIALOPERATORS_H_
#define SEISSOL_SRC_WAVEPROPAGATION_SOLVER_MATERIALOPERATORS_H_

#include "Config.h"
#include "Core/Common/Constants.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/quantities.h"
#include "WavePropagation/Physics/Common.h"
#include "WavePropagation/Physics/Datastructures.h"
#include "WavePropagation/Physics/Quantities.h"
#include "WavePropagation/Solver/Setup.h" // IWYU pragma: export
#include "WavePropagation/Solver/SolverSelector.h"

#include <Eigen/Dense>
#include <array>
#include <complex>
#include <cstddef>

// The operators of a material as the solver of a configuration lays them out, which join the
// physics of the materials (MaterialSetup) to the solvers (SolverSetup). The physics describes a
// material without a solver; the specializations of SolverSetup that this header gathers say how
// each solver lays it out.

namespace seissol::model {

namespace detail {

/// True if the groups have the kinds the codegen laid the tensors out for.
template <std::size_t N, std::size_t M>
constexpr bool kindsMatch(const std::array<QuantityGroup, N>& groups,
                          const std::array<QuantityKind, M>& kinds) {
  if constexpr (N != M) {
    return false;
  } else {
    for (std::size_t i = 0; i < N; ++i) {
      if (groups[i].kind != kinds[i]) {
        return false;
      }
    }
    return true;
  }
}

} // namespace detail

// The codegen and the material headers declare their groups independently --
// the generator cannot read C++, and the material headers are not generated.
// This is what keeps the two from drifting apart: a disagreement becomes a
// compile error rather than a rotation matrix with blocks in the wrong places.
// the argument is a type, which the check takes for an expression in template arguments
// NOLINTBEGIN(bugprone-macro-parentheses)
#define SEISSOL_CHECK_GROUPS(Cfg)                                                                  \
  static_assert(detail::kindsMatch(MaterialOf<Cfg>::RotationGroups<kernels::SolverOf<Cfg>>,        \
                                   generated::Quantities<Cfg>::RotationGroupKinds),                \
                "the material's quantity groups disagree with the generated layout");              \
  static_assert(detail::kindsMatch(MaterialOf<Cfg>::InverseRotationGroups<kernels::SolverOf<Cfg>>, \
                                   generated::Quantities<Cfg>::InverseRotationGroupKinds),         \
                "the material's inverse quantity groups disagree with the generated layout");
// NOLINTEND(bugprone-macro-parentheses)
SEISSOL_FOR_EACH_CONFIG(SEISSOL_CHECK_GROUPS)
#undef SEISSOL_CHECK_GROUPS

/// The operators of a material as the solver of the configuration `Cfg` lays them out.
template <typename Cfg, typename Tmaterial, typename Tmatrix>
void getTransposedCoefficientMatrix(const Tmaterial& material, unsigned dim, Tmatrix& matM) {
  SolverSetup<kernels::SolverOf<Cfg>, Tmaterial>::getTransposedCoefficientMatrix(
      material, dim, matM);
}

template <typename Cfg, typename Tmaterial, typename T>
void getTransposedSourceCoefficientTensor(const Tmaterial& material, T& mE) {
  SolverSetup<kernels::SolverOf<Cfg>, Tmaterial>::getTransposedSourceCoefficientTensor(material,
                                                                                       mE);
}

template <typename Cfg, typename T>
void getPlaneWaveOperator(const T& material,
                          const double n[3],
                          std::complex<double> mdata[T::NumQuantities * T::NumQuantities]) {
  SolverSetup<kernels::SolverOf<Cfg>, T>::getPlaneWaveOperator(material, n, mdata);
}

/// Prepares the data the solver of the configuration `Cfg` keeps per cell.
template <typename Cfg, typename T>
void initializeSpecificLocalData(const T& material,
                                 double timeStepWidth,
                                 typename kernels::SolverOf<Cfg>::LocalData* localData) {
  SolverSetup<kernels::SolverOf<Cfg>, T>::initializeSpecificLocalData(
      material, timeStepWidth, localData);
}

template <typename Cfg, typename T>
void initializeSpecificNeighborData(const T& material,
                                    typename kernels::SolverOf<Cfg>::NeighborData* neighborData) {
  SolverSetup<kernels::SolverOf<Cfg>, T>::initializeSpecificNeighborData(material, neighborData);
}

/**
 * Assembles the face rotation from the material's group declaration.
 *
 * T and Tinv do not always span the same quantities. A solver that keeps the
 * mechanism index in a separate tensor dimension needs the forward rotation to
 * reach one anelastic block, because the flux solver contracts over it, but
 * never applies the inverse there -- so Tinv covers the primary quantities
 * alone and is correspondingly smaller. Each therefore follows its own group
 * list.
 *
 * The inverse is not the transpose for a symmetric second-order tensor, the
 * Voigt weights differ, so each kind supplies a forward and an inverse writer.
 */
template <typename Cfg>
void getFaceRotationMatrix(const CoordinateT& normal,
                           const CoordinateT& tangent1,
                           const CoordinateT& tangent2,
                           typename init::T<Cfg>::view::type& matT,
                           typename init::Tinv<Cfg>::view::type& matTinv) {
  using MaterialT = MaterialOf<Cfg>;
  using SolverT = kernels::SolverOf<Cfg>;
  detail::writeRotationBlocks<false>(
      MaterialT::template RotationGroups<SolverT>, normal, tangent1, tangent2, matT);
  detail::writeRotationBlocks<true>(
      MaterialT::template InverseRotationGroups<SolverT>, normal, tangent1, tangent2, matTinv);
}

template <typename Cfg>
void getFaceRotationMatrix(const Eigen::Vector3d& normal,
                           const Eigen::Vector3d& tangent1,
                           const Eigen::Vector3d& tangent2,
                           typename init::T<Cfg>::view::type& matT,
                           typename init::Tinv<Cfg>::view::type& matTinv) {
  const CoordinateT n = {normal(0), normal(1), normal(2)};
  const CoordinateT s = {tangent1(0), tangent1(1), tangent1(2)};
  const CoordinateT t = {tangent2(0), tangent2(1), tangent2(2)};
  getFaceRotationMatrix<Cfg>(n, s, t, matT, matTinv);
}

} // namespace seissol::model

#endif // SEISSOL_SRC_WAVEPROPAGATION_SOLVER_MATERIALOPERATORS_H_
