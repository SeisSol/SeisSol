// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_EQUATIONS_VISCOACOUSTIC_MODEL_ATTENUATION_H_
#define SEISSOL_SRC_EQUATIONS_VISCOACOUSTIC_MODEL_ATTENUATION_H_

#include <Eigen/Dense>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace seissol::model {
template <std::size_t>
struct ViscoAcousticMaterial;
} // namespace seissol::model

namespace seissol::physics {

template <std::size_t Mechanisms>
void fitAttenuation(seissol::model::ViscoAcousticMaterial<Mechanisms>& vm,
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

    for (size_t i = 0; i < KMax; ++i) {
      for (size_t j = 0; j < Mechanisms; ++j) {
        const double wjsq = w(2 * j) * w(2 * j);
        const double wisq = w(i) * w(i);
        const double norm = wjsq + wisq;
        const double sc1 = w(2 * j) * w(i);
        matAP(i, j) = (sc1 + wjsq / vm.qp) / norm;
      }
    }

    // conversions are to silence NVHPC
    const Eigen::VectorXd qpinv =
        Eigen::VectorXd::Constant(static_cast<std::int64_t>(KMax), 1 / vm.qp);

    auto matAPodc = matAP.completeOrthogonalDecomposition();

    Eigen::VectorXd alpha = matAPodc.solve(qpinv);

    double psi1p = 1;
    double psi2p = 0;
    for (size_t i = 0; i < Mechanisms; ++i) {
      const double w0dw = w0 / w(2 * i);
      const double w0dwsq1 = 1 + w0dw * w0dw;
      psi1p = psi1p - alpha(i) / w0dwsq1;
      psi2p = psi2p + alpha(i) * (w0dw / w0dwsq1);
    }
    const double rp = std::sqrt(psi1p * psi1p + psi2p * psi2p);
    const double varP = vm.lambda * (rp + psi1p) / (2 * rp * rp);

    const double varLambda = varP;
    for (size_t i = 0; i < Mechanisms; ++i) {
      const double t1 = -varP * alpha(i);
      vm.theta[i][0] = t1;
    }

    vm.lambda = varLambda;
  }
}

} // namespace seissol::physics

#endif // SEISSOL_SRC_EQUATIONS_VISCOACOUSTIC_MODEL_ATTENUATION_H_
