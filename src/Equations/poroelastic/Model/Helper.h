// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_HELPER_H_
#define SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_HELPER_H_

#include "Equations/EnergyBase.h"
#include "Equations/elastic/Model/Setup.h"
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

inline AdditionalPoroelasticParameters
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
