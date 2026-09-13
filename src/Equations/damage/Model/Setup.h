// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_SETUP_H_
#define SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_SETUP_H_

#include "Equations/damage/Model/Datastructures.h"
#include "Model/Common.h"

#include <array>
#include <cstddef>
#include <utils/logger.h>

namespace seissol::model {

/// Linearisation of the flux about the cell's mean state.
///
/// This is the operator the Cauchy-Kovalevskaya recursion runs on, not the
/// constitutive law: the law itself is evaluated pointwise by the kernels. The
/// moduli are therefore the effective ones at the cell mean, passed in, and
/// the structure is that of an elastic material embedded in eleven quantities
/// -- the two internal variables carry no flux.
template <typename T>
void getTransposedCoefficientMatrix(
    const DamageMaterial& material, unsigned dim, double lambdaEff, double muEff, T& matM) {
  matM.setZero();

  const double rhoInv = 1.0 / material.rho;

  // The strain equations transport the symmetric velocity gradient. Voigt
  // order is xx, yy, zz, xy, yz, xz, so the shear row of a direction is the
  // one naming it, and it takes a half from each of the two velocities.
  constexpr std::size_t Normal[3] = {0, 1, 2};
  constexpr std::size_t Shear[3][2] = {{3, 5}, {3, 4}, {4, 5}};
  constexpr std::size_t ShearVelocity[3][2] = {{1, 2}, {0, 2}, {1, 0}};

  const std::size_t v = 6 + dim;
  matM(v, Normal[dim]) = -1.0;
  for (std::size_t k = 0; k < 2; ++k) {
    matM(6 + ShearVelocity[dim][k], Shear[dim][k]) = -0.5;
  }

  // The momentum equations transport the traction of a face normal to the
  // direction. Only the isotropic part couples across components.
  for (std::size_t i = 0; i < 3; ++i) {
    matM(i, 6 + dim) = -lambdaEff * rhoInv;
  }
  matM(dim, 6 + dim) -= 2.0 * muEff * rhoInv;
  for (std::size_t k = 0; k < 2; ++k) {
    matM(Shear[dim][k], 6 + ShearVelocity[dim][k]) = -2.0 * muEff * rhoInv;
  }
}

/// What the generic setup asks of this material.
///
/// The linearisation the recursion runs on is the undamaged one: the star
/// matrices are built once per cell, and a linearisation that followed the
/// damage would have to be rebuilt every step. What follows the damage is the
/// constitutive law the kernels evaluate pointwise, and the flux the faces are
/// scaled with -- not this operator.
template <>
struct MaterialSetup<DamageMaterial> : public MaterialSetupDefaults<DamageMaterial> {
  template <typename T>
  static void
      getTransposedCoefficientMatrix(const DamageMaterial& material, unsigned dim, T& matM) {
    seissol::model::getTransposedCoefficientMatrix(
        material, dim, material.lambda0, material.mu0, matM);
  }

  template <typename Tloc, typename Tneigh>
  static void getTransposedGodunovState(const DamageMaterial& /*local*/,
                                        const DamageMaterial& /*neighbor*/,
                                        FaceType /*faceType*/,
                                        Tloc& /*qGodLocal*/,
                                        Tneigh& /*qGodNeighbor*/) {
    // A Godunov state maps the state onto itself, and what this material
    // transports is wider than its state; its faces are scaled from the flux
    // instead, where the pair of matrices is built.
    logError() << "A damaged material has no Godunov state.";
  }

  static DamageMaterial
      getRotatedMaterialCoefficients(const std::array<double, 36>& /*rotationParameters*/,
                                     const DamageMaterial& material) {
    // Isotropic in its undamaged moduli, and the damage is a scalar: nothing
    // to rotate.
    return material;
  }
};

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_SETUP_H_
