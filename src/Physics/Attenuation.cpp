#include <cmath>
#include <cstddef>
#include <array>

#include <Eigen/Dense>

#include "Equations/datastructures.hpp"

namespace seissol::physics {

void fitAttenuation(seissol::model::ViscoElasticMaterial& vm,
                    double freqCentral,
                    double freqRatio) {
#if NUMBER_OF_RELAXATION_MECHANISMS > 0
  constexpr std::size_t nummech = NUMBER_OF_RELAXATION_MECHANISMS;
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
#endif
}

#if defined(USE_ANISOTROPY) && NUMBER_OF_RELAXATION_MECHANISMS > 0
#error                                                                                             \
    "The attenuation still needs to be ported to this case (anisotropy+anelasticity). It could not be done at the time of writing the code because a respective model was not yet established within SeisSol."
#endif

} // namespace seissol::physics
