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
#include <cassert>
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
 * Whether the flux of a face of this type has an operator for the side of the face outside the cell
 * (NeighboringIntegrationData::nAmNm1): the one the corrector applies to the neighbor's trace, or
 * the one a boundary condition applies to the state it builds outside the cell. A free surface
 * builds none, and its Godunov state for that side is a signalling NaN that nothing may read
 * (getTransposedFreeSurfaceGodunovState). Its corrector reads a trace of zeros as the neighbor's
 * all the same, so the operator it applies to it there is zero.
 */
constexpr bool hasExteriorOperator(FaceType faceType) { return faceType != FaceType::FreeSurface; }

/**
 * The transposed coefficient matrix of the face normal a side lifts its flux with: the one of its
 * material, or the one without material where the cell applies the inverse weight of its medium
 * after the lift (see WeightAdjustable).
 */
template <typename MaterialT, typename MatrixT>
void liftedCoefficientMatrix(const MaterialT& material, bool weightAdjusted, MatrixT& matrix) {
  if constexpr (WeightAdjustable) {
    if (weightAdjusted) {
      seissol::model::getTransposedUnweightedCoefficientMatrix(material, 0, matrix);
      return;
    }
  }
  seissol::model::getTransposedCoefficientMatrix(material, 0, matrix);
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
 *
 * A weight-adjusted cell contracts the state with the coefficient matrix without material instead
 * (liftedCoefficientMatrix): the Riemann problem is the same, and the material of the cell comes
 * in at its sample points after the lift. The Godunov state alone has that form; the Rusanov
 * penalty acts on the state and would have to be weighed with the material at the face, so such a
 * cell takes the Godunov flux.
 */
template <typename MaterialT>
void fluxScalarsOfNode(const MaterialT& local,
                       const MaterialT& neighbor,
                       FaceType faceType,
                       parameters::NumericalFlux flux,
                       double fluxScale,
                       std::array<double, FluxCoefficientCount>& plus,
                       std::array<double, FluxCoefficientCount>& minus,
                       bool weightAdjusted = false) {
  assert(!weightAdjusted || flux == parameters::NumericalFlux::Godunov);
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
  liftedCoefficientMatrix(local, weightAdjusted, star);

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

/**
 * The one set of scalars the corrector of a face applies, to the jump of the two traces of the face
 * -- the neighbor's less the cell's own -- out of the two parts of the flux fluxScalarsOfNode gives
 * (see the corrector in the generator). A regular face takes the neighbor's part. A face without a
 * neighbor reads a zero trace as the neighbor's, and takes the cell's own part negated, which then
 * applies the cell's own part to its trace. For a regular face the two agree: the corrector is the
 * strong form wherever a face carries these scalars, and there the cell's own part is the Godunov
 * state less the identity, which is the neighbor's part negated.
 */
template <typename ArrayT>
void correctorScalars(FaceType faceType, const ArrayT& plus, const ArrayT& minus, ArrayT& scalars) {
  static_assert(!NodalFlux || StrongCorrector,
                "the corrector applies one set of scalars to the jump only in the strong form");
  for (std::size_t c = 0; c < scalars.size(); ++c) {
    scalars[c] = faceType == FaceType::Regular ? minus[c] : -plus[c];
  }
}

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_MODEL_CELLFLUX_H_
