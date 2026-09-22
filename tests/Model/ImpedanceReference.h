// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_IMPEDANCEREFERENCE_H_
#define SEISSOL_TESTS_MODEL_IMPEDANCEREFERENCE_H_

#include "Equations/ImpedanceBase.h"
#include "Equations/Setup.h" // IWYU pragma: keep
#include "Model/Common.h"

#include <array>
#include <cstddef>

namespace seissol::unit_test {

/**
 * Y = R_u R_t^-1 from the eigendecomposition of the normal Jacobian: the reference the closed
 * forms of seissol::model::ImpedanceCompute are checked against. Needs MaterialSetup<MaterialT>,
 * i.e. a build which includes the setup of that material.
 *
 * If `lateralStress` is given, the reconstruction R_l R_t^-1 of the remaining stress components is
 * extracted from the very same decomposition -- both express the state jump in the span of the
 * outgoing modes, parametrised by the traction, so the two are consistent by construction.
 */
template <typename MaterialT>
typename seissol::model::ImpedanceCompute<MaterialT>::Matrix admittanceFromEigendecomposition(
    const MaterialT& materialLocal,
    typename seissol::model::ImpedanceCompute<MaterialT>::LateralMatrix* lateralStress = nullptr) {
  using ImpedanceCompute = seissol::model::ImpedanceCompute<MaterialT>;
  using Matrix = typename ImpedanceCompute::Matrix;
  using LateralMatrix = typename ImpedanceCompute::LateralMatrix;

  auto eigenpair = seissol::model::getEigenDecomposition(materialLocal);

  // sigma_ss, sigma_dd, sigma_sd
  const std::array<std::size_t, 3> lateralIndices{1, 2, 4};

  // the eigenvalues are sorted ascendingly by their real part, so the first Dim columns are the
  // modes with a negative wave speed
  const auto matrix = eigenpair.getVectorsAsMatrix();
  Matrix matRT;
  Matrix matRU;
  LateralMatrix matRL;
  for (std::size_t j = 0; j < ImpedanceCompute::Dim; ++j) {
    for (std::size_t i = 0; i < ImpedanceCompute::Dim; ++i) {
      matRT(i, j) = matrix(ImpedanceCompute::TractionIndices[i], j).real();
      matRU(i, j) = matrix(ImpedanceCompute::VelocityIndices[i], j).real();
    }
    for (std::size_t i = 0; i < 3; ++i) {
      matRL(i, j) = matrix(lateralIndices[i], j).real();
    }
  }

  const Matrix matRTInv = matRT.inverse();
  if (lateralStress != nullptr) {
    *lateralStress = matRL * matRTInv;
  }
  return matRU * matRTInv;
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_MODEL_IMPEDANCEREFERENCE_H_
