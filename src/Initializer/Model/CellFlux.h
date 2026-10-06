// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_MODEL_CELLFLUX_H_
#define SEISSOL_SRC_INITIALIZER_MODEL_CELLFLUX_H_

#include "Alignment.h"
#include "Equations/Setup.h" // IWYU pragma: keep
#include "GeneratedCode/coefficients.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Kernels/Precision.h"
#include "Model/Common.h"
#include "Model/OperatorLayout.h"

#include <algorithm>
#include <array>
#include <cstddef>

namespace seissol::initializer {

/**
 * Turns the Godunov state of the local side of a face into the one the corrector of this build
 * applies.
 *
 * Where the operator is one per cell, the corrector is the weak form and the local flux applies
 * the Godunov state as the Riemann problem gives it. Where it varies inside the cell, the volume
 * term is the strong form (see StrongCorrector), and the local flux then subtracts the normal flux
 * of the cell's own trace: the identity comes off the Godunov state. Only the diagonal the state
 * stores is touched -- a solver that folds the relaxation into its quantities keeps the elastic
 * rows alone, and what the relaxation contributes to the normal flux is read through the
 * coefficient matrix the state is contracted with.
 */
template <typename ViewT>
void toCorrectorForm(ViewT& godunovLocal) {
  if constexpr (StrongCorrector) {
    constexpr std::size_t Diagonal =
        std::min(tensor::QgodLocal::Shape[0], tensor::QgodLocal::Shape[1]);
    for (std::size_t i = 0; i < Diagonal; ++i) {
      if (godunovLocal.isInRange(i, i)) {
        godunovLocal(i, i) -= 1;
      }
    }
  }
}

/**
 * The scalars the flux operator of one node is built from, in the coordinates of the face: ten for
 * the Godunov flux of an elastic medium, and one more for the Rusanov penalty on the part of the
 * diagonal the Godunov state never reads, which is zero for the Godunov flux.
 *
 * The matrix form folds the rotation into what a face stores; here it stays in the kernel, because
 * rotated the operator no longer has ten degrees of freedom but fifty-eight. What is stored is the
 * operator as the Riemann problem states it, read at the positions the generated table names.
 *
 * The shape mirrors computeFluxSolverLocal and its neighbour exactly, down to both sides
 * contracting the Godunov state with the coefficient matrix of the *local* material, and the scale
 * the matrix form carries in AplusT riding on the scalars instead.
 *
 * A fault face has no Riemann problem of its own here -- the fault supplies its flux -- so its
 * state is zero. What the local side then still carries is the subtraction of the strong form,
 * see toCorrectorForm, and nothing at all for the weak one.
 */
template <typename MaterialT>
void fluxScalarsOfNode(const MaterialT& local,
                       const MaterialT& neighbor,
                       FaceType faceType,
                       parameters::NumericalFlux flux,
                       double fluxScale,
                       std::array<double, FluxCoefficientCount>& plus,
                       std::array<double, FluxCoefficientCount>& minus) {
  // the Riemann problem is stated over the quantities the Godunov state spans,
  // which is not the count the material declares where a solver carries the
  // relaxation in the same matrix
  constexpr std::size_t N = tensor::QgodLocal::Shape[1];
  constexpr std::size_t Diagonal =
      std::min(tensor::QgodLocal::Shape[0], tensor::QgodLocal::Shape[1]);

  alignas(Alignment) std::array<real, tensor::QgodLocal::size()> godLocalData{};
  alignas(Alignment) std::array<real, tensor::QgodNeighbor::size()> godNeighborData{};
  auto godLocal = init::QgodLocal::view::create(godLocalData.data());
  auto godNeighbor = init::QgodNeighbor::view::create(godNeighborData.data());

  alignas(Alignment) std::array<real, tensor::star::size(0)> starData{};
  auto star = init::star::view<0>::create(starData.data());
  seissol::model::getTransposedCoefficientMatrix(local, 0, star);

  // the Riemann problem, or the central flux the Rusanov form uses instead
  double correction = 0.0;
  if (faceType == FaceType::DynamicRupture) {
    godLocal.setZero();
    godNeighbor.setZero();
  } else if (flux == parameters::NumericalFlux::Rusanov) {
    godLocal.setZero();
    godNeighbor.setZero();
    // the diagonal the Godunov state has: a solver that folds the relaxation
    // into its quantities keeps only the elastic rows of the state
    for (std::size_t i = 0; i < Diagonal; ++i) {
      if (godLocal.isInRange(i, i)) {
        godLocal(i, i) = 0.5;
        godNeighbor(i, i) = 0.5;
      }
    }
    correction = std::max(local.getMaxWaveSpeed(), neighbor.getMaxWaveSpeed()) * 0.5;
  } else {
    seissol::model::getTransposedGodunovState(local, neighbor, faceType, godLocal, godNeighbor);
  }
  toCorrectorForm(godLocal);

  const auto read =
      [&](auto& godunov, double correctionSign, std::array<double, FluxCoefficientCount>& target) {
        for (std::size_t c = 0; c < FluxCoefficientCount; ++c) {
          const auto& source = generated::FluxCoefficientSources[c];
          double value = 0.0;
          for (std::size_t k = 0; k < N; ++k) {
            const double g = godunov.isInRange(source.row, k) ? godunov(source.row, k) : 0.0;
            const double a = star.isInRange(k, source.column) ? star(k, source.column) : 0.0;
            value += g * a;
          }
          // Qcorr is the diagonal the Rusanov form adds, and zero otherwise
          if (source.row == source.column && source.row < Diagonal &&
              godunov.isInRange(source.row, source.row)) {
            value += correctionSign * correction;
          }
          target[c] = fluxScale * value;
        }
      };
  read(godLocal, 1.0, plus);
  read(godNeighbor, -1.0, minus);
}

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_MODEL_CELLFLUX_H_
