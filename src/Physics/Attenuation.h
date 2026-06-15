// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PHYSICS_ATTENUATION_H_
#define SEISSOL_SRC_PHYSICS_ATTENUATION_H_

#include <Eigen/Dense>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace seissol::model {
template <std::size_t>
struct ViscoElasticMaterialParametrized;
} // namespace seissol::model

namespace seissol::physics {

template <std::size_t Mechanisms>
void fitAttenuation(seissol::model::ViscoElasticMaterialParametrized<Mechanisms>& vm,
                    double freqCentral,
                    double freqRatio) {
  if (Mechanisms > 0) {
    constexpr std::size_t KMax =
        2 * Mechanisms - 1; // slight note: if Mechanisms == 0, this does not make any sense

    const double w0 = 2 * M_PI * freqCentral;
    const double wmin = w0 / std::sqrt(freqRatio);

    Eigen::VectorXd w(KMax);

    if (Mechanisms > 1) {
      const double logwmin = std::log(wmin);
      const double logfreqratio = std::log(freqRatio);
      for (std::size_t i = 0; i < KMax; ++i) {
        w(i) = std::exp(logwmin + (i / static_cast<double>(KMax - 1)) * logfreqratio);
      }
    } else {
      w(0) = w0;
    }

    for (size_t i = 0; i < Mechanisms; ++i) {
      vm.omega[i] = w(2 * i);
    }

    Eigen::MatrixXd matAP(KMax, Mechanisms);
    Eigen::MatrixXd matAS(KMax, Mechanisms);

    for (size_t i = 0; i < KMax; ++i) {
      for (size_t j = 0; j < Mechanisms; ++j) {
        const double wjsq = w(2 * j) * w(2 * j);
        const double wisq = w(i) * w(i);
        const double norm = wjsq + wisq;
        const double sc1 = w(2 * j) * w(i);
        matAP(i, j) = (sc1 + wjsq / vm.qp) / norm;
        matAS(i, j) = (sc1 + wjsq / vm.qs) / norm;
      }
    }

    // conversions are to silence NVHPC
    const Eigen::VectorXd qpinv =
        Eigen::VectorXd::Constant(static_cast<std::int64_t>(KMax), 1 / vm.qp);
    const Eigen::VectorXd qsinv =
        Eigen::VectorXd::Constant(static_cast<std::int64_t>(KMax), 1 / vm.qs);

    auto matAPodc = matAP.completeOrthogonalDecomposition();
    auto matASodc = matAS.completeOrthogonalDecomposition();

    Eigen::VectorXd alpha = matAPodc.solve(qpinv);
    Eigen::VectorXd beta = matASodc.solve(qsinv);

    // for anisotropy, use the following:
    //    double var_p = (c11+c22+c33) / 3;
    //    double var_mu = (c44+c55+c66) / 3;
    // and replace it with
    double psi1p = 1;
    double psi2p = 0;
    double psi1s = 1;
    double psi2s = 0;
    for (size_t i = 0; i < Mechanisms; ++i) {
      const double w0dw = w0 / w(2 * i);
      const double w0dwsq1 = 1 + w0dw * w0dw;
      psi1p = psi1p - alpha(i) / w0dwsq1;
      psi2p = psi2p + alpha(i) * (w0dw / w0dwsq1);
      psi1s = psi1s - beta(i) / w0dwsq1;
      psi2s = psi2s + beta(i) * (w0dw / w0dwsq1);
    }
    const double rp = std::sqrt(psi1p * psi1p + psi2p * psi2p);
    const double varP = (vm.lambda + 2 * vm.mu) * (rp + psi1p) / (2 * rp * rp);
    const double rs = std::sqrt(psi1s * psi1s + psi2s * psi2s);
    const double varMu = vm.mu * (rs + psi1s) / (2 * rs * rs);
    // replace end

    const double varLambda = varP - 2 * varMu;
    for (size_t i = 0; i < Mechanisms; ++i) {
      const double t1 = -varP * alpha(i);
      const double t2 = -2.0 * varMu * beta(i);
      vm.theta[i][0] = t1;
      vm.theta[i][1] = t1 - t2;
      vm.theta[i][2] = t2;
    }

    vm.mu = varMu;
    vm.lambda = varLambda;
  }
}

} // namespace seissol::physics

#endif // SEISSOL_SRC_PHYSICS_ATTENUATION_H_
