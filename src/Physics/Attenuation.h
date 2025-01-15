// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PHYSICS_ATTENUATION_H_
#define SEISSOL_SRC_PHYSICS_ATTENUATION_H_

#include "Equations/Datastructures.h"
#include <Eigen/Dense>
#include <cmath>
#include <cstddef>

namespace seissol::physics {

template <std::size_t Mechanisms>
void fitAttenuation(seissol::model::ViscoElasticMaterialParametrized<Mechanisms>& vm,
                    double freqCentral,
                    double freqRatio) {
  if (Mechanisms > 0) {
    constexpr std::size_t nummech = Mechanisms;
    constexpr std::size_t kmax =
        2 * nummech - 1; // slight note: if nummech == 0, this does not make any sense

    const double w0 = 2 * M_PI * freqCentral;
    const double wmin = w0 / std::sqrt(freqRatio);

    Eigen::VectorXd w(kmax);

    if (nummech > 1) {
      const double logwmin = std::log(wmin);
      const double logfreqratio = std::log(freqRatio);
      for (std::size_t i = 0; i < kmax; ++i) {
        w(i) = std::exp(logwmin + (i / static_cast<double>(kmax - 1)) * logfreqratio);
      }
    } else {
      w(0) = w0;
    }

    for (size_t i = 0; i < nummech; ++i) {
      vm.omega[i] = w(2 * i);
    }

    Eigen::MatrixXd AP(kmax, nummech);
    Eigen::MatrixXd AS(kmax, nummech);

    for (size_t i = 0; i < kmax; ++i) {
      for (size_t j = 0; j < nummech; ++j) {
        const double wjsq = w(2 * j) * w(2 * j);
        const double wisq = w(i) * w(i);
        const double norm = wjsq + wisq;
        const double sc1 = w(2 * j) * w(i);
        AP(i, j) = (sc1 + wjsq / vm.Qp) / norm;
        AS(i, j) = (sc1 + wjsq / vm.Qs) / norm;
      }
    }

    Eigen::VectorXd qpinv = Eigen::VectorXd::Constant(kmax, 1 / vm.Qp);
    Eigen::VectorXd qsinv = Eigen::VectorXd::Constant(kmax, 1 / vm.Qs);

    auto APodc = AP.completeOrthogonalDecomposition();
    auto ASodc = AS.completeOrthogonalDecomposition();

    Eigen::VectorXd alpha = APodc.solve(qpinv);
    Eigen::VectorXd beta = ASodc.solve(qsinv);

    // for anisotropy, use the following:
    //    double var_p = (c11+c22+c33) / 3;
    //    double var_mu = (c44+c55+c66) / 3;
    // and replace it with
    double psi1p = 1;
    double psi2p = 0;
    double psi1s = 1;
    double psi2s = 0;
    for (size_t i = 0; i < nummech; ++i) {
      double w0dw = w0 / w(2 * i);
      double w0dwsq1 = 1 + w0dw * w0dw;
      psi1p = psi1p - alpha(i) / w0dwsq1;
      psi2p = psi2p + alpha(i) * (w0dw / w0dwsq1);
      psi1s = psi1s - beta(i) / w0dwsq1;
      psi2s = psi2s + beta(i) * (w0dw / w0dwsq1);
    }
    const double rp = std::sqrt(psi1p * psi1p + psi2p * psi2p);
    const double var_p = (vm.lambda + 2 * vm.mu) * (rp + psi1p) / (2 * rp * rp);
    const double rs = std::sqrt(psi1s * psi1s + psi2s * psi2s);
    const double var_mu = vm.mu * (rs + psi1s) / (2 * rs * rs);
    // replace end

    const double var_lambda = var_p - 2 * var_mu;
    for (size_t i = 0; i < nummech; ++i) {
      const double t1 = -var_p * alpha(i);
      const double t2 = -2.0 * var_mu * beta(i);
      vm.theta[i][0] = t1;
      vm.theta[i][1] = t1 - t2;
      vm.theta[i][2] = t2;
    }

    vm.mu = var_mu;
    vm.lambda = var_lambda;
  }
}

} // namespace seissol::physics

#endif // SEISSOL_SRC_PHYSICS_ATTENUATION_H_
