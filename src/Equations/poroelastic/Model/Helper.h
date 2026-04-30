// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_HELPER_H_
#define SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_HELPER_H_

#include "Equations/elastic/Model/ElasticSetup.h"
#include "Equations/poroelastic/Model/Datastructures.h"
#include "GeneratedCode/init.h"
#include "Kernels/Common.h"
#include "Model/Common.h"
#include "Numerical/Eigenvalues.h"
#include "Numerical/Transformation.h"

#include <Eigen/Dense>
#include <cassert>
#include <yateto.h>

namespace seissol::init {
class Z;
class Zinv;
} // namespace seissol::init

namespace seissol::model {

template <typename Tview>
static void calcZinv(yateto::DenseTensorView<2, real, unsigned>& zInv,
                     Tview& sourceMatrix,
                     size_t quantity,
                     double timeStepWidth) {
  using Matrix = Eigen::Matrix<real, ConvergenceOrder, ConvergenceOrder>;
  using Vector = Eigen::Matrix<real, ConvergenceOrder, 1>;

  Matrix matZ{init::Z::Values};
  // sourceMatrix[i,i] = 0 for i < 10
  // This is specific to poroelasticity, so change this for another equation
  // We need this check, because otherwise the lookup sourceMatrix(quantity, quantity) fails
  if (quantity >= 10) {
    matZ -= timeStepWidth * sourceMatrix(quantity, quantity) * Matrix::Identity();
  }

  auto solver = matZ.colPivHouseholderQr();
  for (std::size_t col = 0; col < ConvergenceOrder; col++) {
    Vector rhs = Vector::Zero();
    rhs(col) = 1.0;
    auto zInvCol = solver.solve(rhs);
    for (std::size_t row = 0; row < ConvergenceOrder; row++) {
      // save as transposed
      zInv(col, row) = zInvCol(row);
    }
  }
}

// constexpr for loop since we need to instatiate the view templates
template <size_t Istart, size_t Iend, typename Tview>
struct ZInvInitializer {
  ZInvInitializer(
      real zInvData[PoroElasticMaterial::NumQuantities][ConvergenceOrder * ConvergenceOrder],
      Tview& sourceMatrix,
      real timeStepWidth) {
    auto zInv = init::Zinv::view<Istart>::create(zInvData[Istart]);
    calcZinv(zInv, sourceMatrix, Istart, timeStepWidth);
    if constexpr (Istart < Iend - 1) {
      ZInvInitializer<Istart + 1, Iend, Tview>(zInvData, sourceMatrix, timeStepWidth);
    }
  };
};

struct AdditionalPoroelasticParameters {
  Eigen::Matrix<double, 6, 1> alpha;
  // NOLINTNEXTLINE
  double KBar;
  // NOLINTNEXTLINE
  double M;
  double m;
  Eigen::Matrix<double, 6, 6> cBar;
  double rhoBar;
  double rho1;
  double rho2;
  double beta1;
  double beta2;
};

static AdditionalPoroelasticParameters
    getAdditionalParameters(const PoroElasticMaterial& material) {
  Eigen::Matrix<double, 6, 1> alpha;
  alpha << 1 - (3 * material.lambda + 2 * material.mu) / (3 * material.bulkSolid),
      1 - (3 * material.lambda + 2 * material.mu) / (3 * material.bulkSolid),
      1 - (3 * material.lambda + 2 * material.mu) / (3 * material.bulkSolid), -0.0, -0.0, -0.0;

  Eigen::Matrix<double, 6, 6> c;
  c << material.lambda + 2 * material.mu, material.lambda, material.lambda, 0, 0, 0,
      material.lambda, material.lambda + 2 * material.mu, material.lambda, 0, 0, 0, material.lambda,
      material.lambda, material.lambda + 2 * material.mu, 0, 0, 0, 0, 0, 0, material.mu, 0, 0, 0, 0,
      0, 0, material.mu, 0, 0, 0, 0, 0, 0, material.mu;

  const double cKBar = material.lambda + 2 * material.mu / 3;
  const double cM =
      material.bulkSolid / (1 - material.porosity - cKBar / material.bulkSolid +
                            material.porosity * material.bulkSolid / material.bulkFluid);
  const double m = material.rhoFluid * material.tortuosity / material.porosity;

  Eigen::Matrix<double, 6, 6> cBar = c + cM * alpha * alpha.transpose();

  const double rhoBar =
      (1 - material.porosity) * material.rho + material.porosity * material.rhoFluid;
  const double rho1 = rhoBar - material.rhoFluid * material.rhoFluid / m;
  const double rho2 = material.rhoFluid - m * rhoBar / material.rhoFluid;
  const double beta1 = material.rhoFluid / m;
  const double beta2 = rhoBar / material.rhoFluid;

  return {alpha, cKBar, cM, m, cBar, rhoBar, rho1, rho2, beta1, beta2};
}

} // namespace seissol::model
#endif // SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_HELPER_H_
