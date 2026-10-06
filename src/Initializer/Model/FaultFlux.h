// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_MODEL_FAULTFLUX_H_
#define SEISSOL_SRC_INITIALIZER_MODEL_FAULTFLUX_H_

#include "Alignment.h"
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Typedefs.h"
#include "Equations/Setup.h" // IWYU pragma: keep
#include "GeneratedCode/coefficients.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Precision.h"
#include "Model/Common.h"
#include "Model/CommonDatastructures.h"
#include "Model/OperatorLayout.h"
#include "Solver/MultipleSimulations.h"

#include <algorithm>
#include <array>
#include <cstddef>

namespace seissol::initializer {

/**
 * The lift of one side of a fault face, where the face carries it as one matrix (see
 * dr::FaultFluxLayout); `target` takes tensor::fluxSolver::size() reals.
 *
 * The lift applies the coefficient matrix of the fault normal to the imposed state, which is given
 * in the coordinates of the face. So it is the star of the first direction of the material seen in
 * those coordinates -- the rotated material, which the impedance matrices of an anisotropic or a
 * poroelastic face and the flux of a regular face are formed from as well; rotating an isotropic
 * material changes nothing. The scale of the side and the rotation back to global coordinates are
 * folded into the matrix.
 */
template <typename MaterialT>
void setMatrixFaultFlux(real* target,
                        const real* rotation,
                        double fluxScale,
                        const MaterialT& material,
                        const std::array<double, 36>& bond) {
  alignas(Alignment) std::array<real, tensor::star::size(0)> star{};
  auto viewStar = init::star::view<0>::create(star.data());
  seissol::model::getTransposedCoefficientMatrix(
      seissol::model::getRotatedMaterialCoefficients(bond, material), 0, viewStar);

  dynamicRupture::kernel::rotateFluxMatrix krnl;
  krnl.T = rotation;
  krnl.fluxSolver = target;
  krnl.fluxScaleDR = fluxScale;
  krnl.star(0) = star.data();
  krnl.execute();
}

/**
 * The lift of one side of a fault face, where the face carries it per point (see
 * dr::FaultFluxLayout).
 *
 * The lift applies the coefficient matrix of the fault normal to the imposed state, which is given
 * in the coordinates of the face. So at every point it is the star of the first direction of the
 * material there, turned into those coordinates -- as setMatrixFaultFlux takes it, and as the
 * impedance matrices of that point take it where the material has them. The rotation back is
 * stored once, and the scale of the side rides on the scalars, just as the matrix form folds both
 * into one matrix.
 *
 * A scalar that the solver does not read off the material -- a relaxation frequency, or the unit
 * weight of a relaxation block -- is not among what the fault evaluates at its points, and it is
 * one number for the whole domain anyway, so it is taken from the cell.
 *
 * `atPoints` holds the material at the points of the face in the order the impedances read it,
 * with the fused simulations innermost.
 */
template <typename MaterialT, std::size_t Points>
void setPointwiseFaultFlux(real* target,
                           const real* rotation,
                           const std::array<double, dr::misc::NumBoundaryGaussPoints>& fluxScale,
                           const std::array<MaterialT, Points>& atPoints,
                           const MaterialT& cellMaterial,
                           const std::array<double, 36>& bond) {
  using Setup = seissol::model::SolverSetup<typename MaterialT::Solver, MaterialT>;
  // a build that keeps one matrix per side never calls this, but still compiles it
  static_assert(!NodalFaultFlux ||
                    Points >= dr::misc::NumBoundaryGaussPoints * multisim::NumSimulations,
                "The fault lift reads more points than the material is given at.");

  std::fill_n(target, dr::FaultFluxLayout::Size, static_cast<real>(0));
  std::copy_n(
      rotation, dr::FaultFluxLayout::RotationSize, target + dr::FaultFluxLayout::RotationOffset);

  const auto ofCell = seissol::model::getStarCoefficients(cellMaterial);
  for (std::size_t point = 0; point < dr::misc::NumBoundaryGaussPoints; ++point) {
    // the material is one field, so every fused simulation sees the same at a point
    const auto ofPoint =
        seissol::model::getStarCoefficients(seissol::model::getRotatedMaterialCoefficients(
            bond, atPoints[point * multisim::NumSimulations]));
    for (std::size_t coefficient = 0; coefficient < FaultFluxCoefficientCount; ++coefficient) {
      const auto index = generated::FaultFluxCoefficientIndices[coefficient];
      const double value =
          Setup::CoefficientOrigins[index] == seissol::model::CoefficientOrigin::Material
              ? ofPoint[index]
              : ofCell[index];
      target[dr::FaultFluxLayout::coefficientOffset(coefficient) + point] =
          static_cast<real>(fluxScale[point] * value);
    }
  }
}

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_MODEL_FAULTFLUX_H_
