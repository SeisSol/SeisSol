// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_SETUP_H_
#define SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_SETUP_H_

#include "Equations/elastic/Model/Setup.h"
#include "Equations/poroelastic/Model/Datastructures.h"
#include "Equations/poroelastic/Model/Helper.h"
#include "GeneratedCode/init.h"
#include "Kernels/Common.h"
#include "Model/Common.h"
#include "Numerical/Eigenvalues.h"
#include "Numerical/Transformation.h"

#include <Eigen/Dense>
#include <cassert>
#include <yateto.h>

namespace seissol::model {

#ifdef SEISSOL_KERNELS_STP

template <>
struct MaterialSetup<PoroElasticMaterial> : public MaterialSetupDefaults<PoroElasticMaterial> {
  template <typename T>
  static void setToZero(T& matM) {
    matM.setZero();
  }

  template <typename T>
  static void
      getTransposedCoefficientMatrix(const PoroElasticMaterial& material, unsigned dim, T& matM) {
    setToZero<T>(matM);
    const AdditionalPoroelasticParameters params = getAdditionalParameters(material);
    switch (dim) {
    case 0:
      matM(0, 6) = -1 / params.rho1;
      matM(0, 10) = -1 / params.rho2;
      matM(3, 7) = -1 / params.rho1;
      matM(3, 11) = -1 / params.rho2;
      matM(5, 8) = -1 / params.rho1;
      matM(5, 12) = -1 / params.rho2;

      matM(6, 0) = -params.cBar(0, 0);
      matM(6, 1) = -params.cBar(1, 0);
      matM(6, 2) = -params.cBar(2, 0);
      matM(6, 3) = -params.cBar(5, 0);
      matM(6, 4) = -params.cBar(3, 0);
      matM(6, 5) = -params.cBar(4, 0);
      matM(6, 9) = params.M * params.alpha(0);

      matM(7, 0) = -params.cBar(0, 5);
      matM(7, 1) = -params.cBar(1, 5);
      matM(7, 2) = -params.cBar(2, 5);
      matM(7, 3) = -params.cBar(5, 5);
      matM(7, 4) = -params.cBar(3, 5);
      matM(7, 5) = -params.cBar(4, 5);
      matM(7, 9) = params.M * params.alpha(5);

      matM(8, 0) = -params.cBar(0, 4);
      matM(8, 1) = -params.cBar(1, 4);
      matM(8, 2) = -params.cBar(2, 4);
      matM(8, 3) = -params.cBar(5, 4);
      matM(8, 4) = -params.cBar(3, 4);
      matM(8, 5) = -params.cBar(4, 4);
      matM(8, 9) = params.M * params.alpha(4);

      matM(9, 6) = -params.beta1 / params.rho1;
      matM(9, 10) = -params.beta2 / params.rho2;

      matM(10, 0) = -params.M * params.alpha(0);
      matM(10, 1) = -params.M * params.alpha(1);
      matM(10, 2) = -params.M * params.alpha(2);
      matM(10, 3) = -params.M * params.alpha(5);
      matM(10, 4) = -params.M * params.alpha(3);
      matM(10, 5) = -params.M * params.alpha(4);
      matM(10, 9) = params.M;
      break;
    case 1:
      matM(1, 7) = -1 / params.rho1;
      matM(1, 11) = -1 / params.rho2;
      matM(3, 6) = -1 / params.rho1;
      matM(3, 10) = -1 / params.rho2;
      matM(4, 8) = -1 / params.rho1;
      matM(4, 12) = -1 / params.rho2;

      matM(6, 0) = -params.cBar(0, 5);
      matM(6, 1) = -params.cBar(1, 5);
      matM(6, 2) = -params.cBar(2, 5);
      matM(6, 3) = -params.cBar(5, 5);
      matM(6, 4) = -params.cBar(3, 5);
      matM(6, 5) = -params.cBar(4, 5);
      matM(6, 9) = params.M * params.alpha(5);

      matM(7, 0) = -params.cBar(0, 1);
      matM(7, 1) = -params.cBar(1, 1);
      matM(7, 2) = -params.cBar(2, 1);
      matM(7, 3) = -params.cBar(5, 1);
      matM(7, 4) = -params.cBar(3, 1);
      matM(7, 5) = -params.cBar(4, 1);
      matM(7, 9) = params.M * params.alpha(1);

      matM(8, 0) = -params.cBar(0, 3);
      matM(8, 1) = -params.cBar(1, 3);
      matM(8, 2) = -params.cBar(2, 3);
      matM(8, 3) = -params.cBar(5, 3);
      matM(8, 4) = -params.cBar(3, 3);
      matM(8, 5) = -params.cBar(4, 3);
      matM(8, 9) = params.M * params.alpha(3);

      matM(9, 7) = -params.beta1 / params.rho1;
      matM(9, 11) = -params.beta2 / params.rho2;

      matM(11, 0) = -params.M * params.alpha(0);
      matM(11, 1) = -params.M * params.alpha(1);
      matM(11, 2) = -params.M * params.alpha(2);
      matM(11, 3) = -params.M * params.alpha(5);
      matM(11, 4) = -params.M * params.alpha(3);
      matM(11, 5) = -params.M * params.alpha(4);
      matM(11, 9) = params.M;
      break;
    case 2:
      matM(2, 8) = -1 / params.rho1;
      matM(2, 12) = -1 / params.rho2;
      matM(4, 7) = -1 / params.rho1;
      matM(4, 11) = -1 / params.rho2;
      matM(5, 6) = -1 / params.rho1;
      matM(5, 10) = -1 / params.rho2;

      matM(6, 0) = -params.cBar(0, 4);
      matM(6, 1) = -params.cBar(1, 4);
      matM(6, 2) = -params.cBar(2, 4);
      matM(6, 3) = -params.cBar(5, 4);
      matM(6, 4) = -params.cBar(3, 4);
      matM(6, 5) = -params.cBar(4, 4);
      matM(6, 9) = params.M * params.alpha(4);

      matM(7, 0) = -params.cBar(0, 3);
      matM(7, 1) = -params.cBar(1, 3);
      matM(7, 2) = -params.cBar(2, 3);
      matM(7, 3) = -params.cBar(5, 3);
      matM(7, 4) = -params.cBar(3, 3);
      matM(7, 5) = -params.cBar(4, 3);
      matM(7, 9) = params.M * params.alpha(3);

      matM(8, 0) = -params.cBar(0, 2);
      matM(8, 1) = -params.cBar(1, 2);
      matM(8, 2) = -params.cBar(2, 2);
      matM(8, 3) = -params.cBar(5, 2);
      matM(8, 4) = -params.cBar(3, 2);
      matM(8, 5) = -params.cBar(4, 2);
      matM(8, 9) = params.M * params.alpha(2);

      matM(9, 8) = -params.beta1 / params.rho1;
      matM(9, 12) = -params.beta2 / params.rho2;

      matM(12, 0) = -params.M * params.alpha(0);
      matM(12, 1) = -params.M * params.alpha(1);
      matM(12, 2) = -params.M * params.alpha(2);
      matM(12, 3) = -params.M * params.alpha(5);
      matM(12, 4) = -params.M * params.alpha(3);
      matM(12, 5) = -params.M * params.alpha(4);
      matM(12, 9) = params.M;
      break;

    default:
      logError() << "Cannot create transposed coefficient matrix for dimension " << dim
                 << ", has to be either 0, 1 or 2.";
    }
  }

  template <typename T>
  static void getTransposedSourceCoefficientTensor(const PoroElasticMaterial& material, T& matE) {
    const AdditionalPoroelasticParameters params = getAdditionalParameters(material);
    const double e1 = params.beta1 * material.viscosity / (params.rho1 * material.permeability);
    const double e2 = params.beta2 * material.viscosity / (params.rho2 * material.permeability);

    matE.setZero();
    matE(10, 6) = e1;
    matE(11, 7) = e1;
    matE(12, 8) = e1;

    matE(10, 10) = e2;
    matE(11, 11) = e2;
    matE(12, 12) = e2;
  }

  template <typename Tloc, typename Tneigh>
  static void getTransposedGodunovState(const PoroElasticMaterial& local,
                                        const PoroElasticMaterial& neighbor,
                                        FaceType faceType,
                                        Tloc& qGodLocal,
                                        Tneigh& qGodNeighbor) {
    // Will be used to check, whether numbers are (numerically) zero
    constexpr auto ZeroThreshold = 1e-7;
    using CMatrix = Eigen::Matrix<std::complex<double>,
                                  PoroElasticMaterial::NumQuantities,
                                  PoroElasticMaterial::NumQuantities>;
    using Matrix = Eigen::
        Matrix<double, PoroElasticMaterial::NumQuantities, PoroElasticMaterial::NumQuantities>;
    using CVector = Eigen::Matrix<std::complex<double>, PoroElasticMaterial::NumQuantities, 1>;

    auto splitEigenDecomposition = [](const PoroElasticMaterial& material) {
      auto eigenpair = getEigenDecomposition(material, ZeroThreshold);
      return std::pair<CVector, CMatrix>{eigenpair.getValuesAsVector(),
                                         eigenpair.getVectorsAsMatrix()};
    };

    auto [localEigenvalues, localEigenvectors] = splitEigenDecomposition(local);
    auto [neighborEigenvalues, neighborEigenvectors] = splitEigenDecomposition(neighbor);

    CMatrix chiMinus = CMatrix::Zero();
    CMatrix chiPlus = CMatrix::Zero();
    for (int i = 0; i < 13; i++) {
      if (localEigenvalues(i).real() < -ZeroThreshold) {
        chiMinus(i, i) = 1.0;
      }
      if (localEigenvalues(i).real() > ZeroThreshold) {
        chiPlus(i, i) = 1.0;
      }
    }

    // matR == eigenvector matrix
    CMatrix matR = localEigenvectors * chiMinus + neighborEigenvectors * chiPlus;
    // set null space eigenvectors manually
    matR(1, 4) = 1.0;
    matR(2, 5) = 1.0;
    matR(12, 6) = 1.0;
    matR(11, 7) = 1.0;
    matR(4, 8) = 1.0;
    if (faceType == FaceType::FreeSurface) {
      Matrix realR = matR.real();
      getTransposedFreeSurfaceGodunovState<PoroElasticMaterial>(
          MaterialType::Poroelastic, qGodLocal, qGodNeighbor, realR);
    } else {
      // Only the outgoing (negative eigenvalue) projector is computed; qGodLocal is its complement.
      //
      // Note that chiMinus and chiPlus are NOT complementary here: the five zero eigenvalues are in
      // neither of them, so forming both projectors separately leaves the null-space modes out of
      // both Godunov matrices. Deriving qGodLocal as I - qGodNeighbor assigns them to the local
      // subsystem, exactly as the elastic and anisotropic setups already do. This does not change
      // the flux at all -- the affected columns (1, 2, 4, 11, 12) multiply the structurally zero
      // rows of the star matrix -- but it restores qGodLocal + qGodNeighbor == I, which is what
      // makes the invariant testable and the flux solver assembly consistent across equation sets.
      const auto matRT = matR.transpose();

      // Deliberately partialPivLu and not a rank-revealing factorization; cf. ElasticSetup.h
      const auto matRlu = matRT.partialPivLu();
      const auto godunovMinus = matRlu.solve(chiMinus * matRT).eval();

      for (unsigned i = 0; i < qGodLocal.shape(0); ++i) {
        for (unsigned j = 0; j < qGodLocal.shape(1); ++j) {
          const double identity = (i == j) ? 1.0 : 0.0;
          qGodLocal(i, j) = identity - godunovMinus(i, j).real();
          qGodNeighbor(i, j) = godunovMinus(i, j).real();
          assert(std::abs(godunovMinus(i, j).imag()) < ZeroThreshold);
        }
      }
    }
  }
};

#endif

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_SETUP_H_
