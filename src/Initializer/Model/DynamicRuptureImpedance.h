// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_MODEL_DYNAMICRUPTUREIMPEDANCE_H_
#define SEISSOL_SRC_INITIALIZER_MODEL_DYNAMICRUPTUREIMPEDANCE_H_

#include "Equations/Datastructures.h" // IWYU pragma: keep
#include "Equations/Setup.h"          // IWYU pragma: keep
#include "GeneratedCode/tensor.h"
#include "Model/Common.h"
#include "Model/CommonDatastructures.h"
#include "Numerical/Eigenvalues.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <cstddef>
#include <optional>
#include <sstream>
#include <string>

namespace seissol::initializer::model {

/// Number of coupled traction/velocity components at a fault: 3 for
/// (visco)elastic and anisotropic, 4 for poroelastic (fluid pressure/flux).
constexpr std::size_t DrImpedanceDim = seissol::tensor::Zminus::Shape[0];
using DrMatrix = Eigen::Matrix<double, DrImpedanceDim, DrImpedanceDim>;
/// Maps the interface traction to the three stress components which do not take part in the
/// fault-normal Riemann problem (sigma_ss, sigma_dd, sigma_sd).
using DrLateralMatrix = Eigen::Matrix<double, 3, DrImpedanceDim>;

/**
 * All impedance-derived quantities of a single fault face, in fault-local
 * coordinates (n, s, d [, fluid]).
 *
 * NOTE ON NAMING: what the `Zplus`/`Zminus` tensors and
 * `dr::ImpedanceMatrices::impedance{,Neig}` store is *not* an impedance but an
 * admittance Y, mapping traction to velocity. The quantity carrying the
 * dimension of an impedance is `eta`.
 *
 *   eta = (Y+ + Y-)^-1
 *   b+  = eta * Y+ ,  b- = eta * Y-   (traction averaging, b+ + b- = I)
 *
 * Isotropic elastic limit:
 *   Y+  = diag(1/Zp, 1/Zs, 1/Zs)
 *   eta = diag(etaP, etaS, etaS)
 *   b+  = diag(etaP/Zp, etaS/Zs, etaS/Zs)
 * i.e. exactly the values the "fast" branch writes into tractionPlusMatrix.
 */
struct FaultImpedance {
  DrMatrix admittancePlus;
  DrMatrix admittanceMinus;
  DrMatrix eta;
  DrMatrix bPlus;
  DrMatrix bMinus;
  /// only the plus side, matching the convention of ReceiverOutput::computeLocalStresses
  DrLateralMatrix lateralStressPlus;
};

namespace impedance_detail {

/**
 * Christoffel matrix for propagation along the local x axis, i.e. along the
 * fault normal once the material has been rotated with the Bond matrix:
 *
 *   Gamma_ik = C_{i1k1}
 *
 * In Voigt notation this is the (11, 16, 15 / 16, 66, 56 / 15, 56, 55) block --
 * the same coefficients that appear in rows 6..8 of
 * getTransposedCoefficientMatrix(material, 0, .).
 */
inline Eigen::Matrix3d christoffelMatrix(const seissol::model::AnisotropicMaterial& materialLocal) {
  Eigen::Matrix3d gamma;
  gamma(0, 0) = materialLocal.c11;
  gamma(0, 1) = materialLocal.c16;
  gamma(0, 2) = materialLocal.c15;
  gamma(1, 1) = materialLocal.c66;
  gamma(1, 2) = materialLocal.c56;
  gamma(2, 2) = materialLocal.c55;
  gamma(1, 0) = gamma(0, 1);
  gamma(2, 0) = gamma(0, 2);
  gamma(2, 1) = gamma(1, 2);
  return gamma;
}

/**
 * Reconstructs the stress components which do not take part in the fault-normal Riemann problem
 * from the traction difference across the fault.
 *
 * A plane wave travelling along the local x axis only has the strain rates
 * (eps_1, eps_6, eps_5) = (du/dx, dv/dx, dw/dx), so in Voigt notation
 *
 *   [d sigma_1; d sigma_6; d sigma_5] = Gamma * [d eps_1; d eps_6; d eps_5]
 *   [d sigma_2; d sigma_3; d sigma_4] = C[{2,3,4},{1,6,5}] * [d eps_1; d eps_6; d eps_5]
 *
 * with the very same Gamma as above. Eliminating the strain rates gives
 *
 *   [d sigma_ss; d sigma_dd; d sigma_sd] = C[{2,3,4},{1,6,5}] * Gamma^-1 * [traction difference].
 *
 * For an isotropic material this reduces to d sigma_ss = d sigma_dd = d sigma_nn * (1 - 2
 * (cs/cp)^2) and d sigma_sd = 0, which is what the fault output used before.
 */
inline Eigen::Matrix3d
    lateralStressFromChristoffel(const seissol::model::AnisotropicMaterial& materialLocal) {
  Eigen::Matrix3d rest;
  // rows: Voigt 2, 3, 4 (sigma_ss, sigma_dd, sigma_sd); columns: Voigt 1, 6, 5
  rest(0, 0) = materialLocal.c12;
  rest(0, 1) = materialLocal.c26;
  rest(0, 2) = materialLocal.c25;
  rest(1, 0) = materialLocal.c13;
  rest(1, 1) = materialLocal.c36;
  rest(1, 2) = materialLocal.c35;
  rest(2, 0) = materialLocal.c14;
  rest(2, 1) = materialLocal.c46;
  rest(2, 2) = materialLocal.c45;
  return rest * christoffelMatrix(materialLocal).inverse();
}

/**
 * Admittance of an elastic/anisotropic half space in closed form.
 *
 *   Z = sum_alpha rho * v_alpha * p_alpha (x) p_alpha = (rho * Gamma)^(1/2)
 *   Y = Z^-1
 *
 * because Gamma p_alpha = rho v_alpha^2 p_alpha with an orthonormal p_alpha.
 * Gamma is symmetric positive definite, so a self-adjoint eigensolver is used
 * and the result is symmetric positive definite by construction -- also when
 * qS1 and qS2 are degenerate, which is exactly where the general 9x9 complex
 * eigendecomposition loses accuracy.
 */
inline Eigen::Matrix3d admittanceFromChristoffel(const Eigen::Matrix3d& gamma, double rho) {
  const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(rho * gamma);
  const Eigen::Vector3d impedanceEigenvalues = solver.eigenvalues().cwiseSqrt();
  return solver.eigenvectors() * impedanceEigenvalues.cwiseInverse().asDiagonal() *
         solver.eigenvectors().transpose();
}

/**
 * Y = R_u R_t^-1 from the eigendecomposition of the normal Jacobian, kept for poroelastic where no
 * closed form is available here.
 *
 * If `lateralStress` is given, the reconstruction R_l R_t^-1 of the remaining stress components is
 * extracted from the very same decomposition -- both express the state jump in the span of the
 * outgoing modes, parametrised by the traction, so the two are consistent by construction.
 */
template <typename MaterialT>
DrMatrix admittanceFromEigendecomposition(const MaterialT& materialLocal,
                                          DrLateralMatrix* lateralStress = nullptr) {
  auto eigenpair = seissol::model::getEigenDecomposition(materialLocal);

  std::vector<std::size_t> tractionIndices;
  std::vector<std::size_t> velocityIndices;
  std::vector<std::size_t> columnIndices;
  if constexpr (MaterialT::Type == seissol::model::MaterialType::Poroelastic) {
    tractionIndices = {0, 3, 5, 9};
    velocityIndices = {6, 7, 8, 10};
    columnIndices = {0, 1, 2, 3};
  } else {
    tractionIndices = {0, 3, 5};
    velocityIndices = {6, 7, 8};
    columnIndices = {0, 1, 2};
  }

  // sigma_ss, sigma_dd, sigma_sd
  const std::array<std::size_t, 3> lateralIndices{1, 2, 4};

  const auto matrix = eigenpair.getVectorsAsMatrix();
  DrMatrix matRT;
  DrMatrix matRU;
  DrLateralMatrix matRL;
  for (std::size_t j = 0; j < DrImpedanceDim; ++j) {
    for (std::size_t i = 0; i < DrImpedanceDim; ++i) {
      matRT(i, j) = matrix(tractionIndices[i], columnIndices[j]).real();
      matRU(i, j) = matrix(velocityIndices[i], columnIndices[j]).real();
    }
    for (std::size_t i = 0; i < 3; ++i) {
      matRL(i, j) = matrix(lateralIndices[i], columnIndices[j]).real();
    }
  }

  const DrMatrix matRTInv = matRT.inverse();
  if (lateralStress != nullptr) {
    *lateralStress = matRL * matRTInv;
  }
  return matRU * matRTInv;
}

} // namespace impedance_detail

/// Admittance of one side. `materialLocal` must already be rotated into the
/// fault-local frame (Bond matrix).
template <typename MaterialT>
DrMatrix computeAdmittance(const MaterialT& materialLocal,
                           DrLateralMatrix* lateralStress = nullptr) {
  if constexpr (MaterialT::Type == seissol::model::MaterialType::Anisotropic) {
    if (lateralStress != nullptr) {
      *lateralStress = impedance_detail::lateralStressFromChristoffel(materialLocal);
    }
    return impedance_detail::admittanceFromChristoffel(
        impedance_detail::christoffelMatrix(materialLocal), materialLocal.rho);
  } else {
    return impedance_detail::admittanceFromEigendecomposition(materialLocal, lateralStress);
  }
}

inline FaultImpedance assembleFaultImpedance(const DrMatrix& admittancePlus,
                                             const DrMatrix& admittanceMinus) {
  FaultImpedance impedance;
  impedance.admittancePlus = admittancePlus;
  impedance.admittanceMinus = admittanceMinus;
  impedance.eta = (admittancePlus + admittanceMinus).inverse();
  impedance.bPlus = impedance.eta * admittancePlus;
  impedance.bMinus = impedance.eta * admittanceMinus;
  return impedance;
}

/// Convenience: both sides at once, materials given in the fault-local frame.
template <typename MaterialT>
FaultImpedance computeFaultImpedance(const MaterialT& plusLocal, const MaterialT& minusLocal) {
  DrLateralMatrix lateralStressPlus = DrLateralMatrix::Zero();
  const auto admittancePlus = computeAdmittance(plusLocal, &lateralStressPlus);
  auto impedance = assembleFaultImpedance(admittancePlus, computeAdmittance(minusLocal));
  impedance.lateralStressPlus = lateralStressPlus;
  return impedance;
}

/**
 * Checks the invariants every fault impedance has to satisfy. Returns
 * std::nullopt on success, otherwise a human readable description of the first
 * violation.
 *
 * `expectSelfAdjoint` should be false for poroelastic, where the 4x4 block is
 * not built from an energy-conjugate pair.
 */
inline std::optional<std::string>
    checkFaultImpedance(const FaultImpedance& impedance,
                        bool expectSelfAdjoint = seissol::model::MaterialT::Type !=
                                                 seissol::model::MaterialType::Poroelastic,
                        double tolerance = 1e-9) {
  const auto relIdentityError = [](const DrMatrix& matrix) {
    return (matrix - DrMatrix::Identity()).cwiseAbs().maxCoeff();
  };
  const auto relSymmetryError = [](const DrMatrix& matrix) {
    const double scale = matrix.cwiseAbs().maxCoeff();
    return (scale == 0.0) ? 0.0 : (matrix - matrix.transpose()).cwiseAbs().maxCoeff() / scale;
  };
  const auto smallestEigenvalue = [](const DrMatrix& matrix) {
    const DrMatrix symmetric = 0.5 * (matrix + matrix.transpose());
    const Eigen::SelfAdjointEigenSolver<DrMatrix> solver(symmetric, Eigen::EigenvaluesOnly);
    return solver.eigenvalues().minCoeff();
  };

  std::ostringstream message;

  if (!impedance.admittancePlus.allFinite() || !impedance.admittanceMinus.allFinite() ||
      !impedance.eta.allFinite() || !impedance.bPlus.allFinite() || !impedance.bMinus.allFinite()) {
    message << "non-finite entry in the fault impedance";
    return message.str();
  }

  // catches a badly conditioned inverse, e.g. from a degenerate eigenbasis
  const double etaError =
      relIdentityError(impedance.eta * (impedance.admittancePlus + impedance.admittanceMinus));
  if (etaError > tolerance) {
    message << "eta * (Y+ + Y-) != I, max error " << etaError;
    return message.str();
  }

  // catches an eta*Z / eta*Y mixup and any inconsistent transpose convention
  const double bError = relIdentityError(impedance.bPlus + impedance.bMinus);
  if (bError > tolerance) {
    message << "b+ + b- != I, max error " << bError;
    return message.str();
  }

  if (expectSelfAdjoint) {
    const double symmetryError = std::max({relSymmetryError(impedance.admittancePlus),
                                           relSymmetryError(impedance.admittanceMinus),
                                           relSymmetryError(impedance.eta)});
    if (symmetryError > tolerance) {
      message << "fault impedance is not self-adjoint, max relative error " << symmetryError;
      return message.str();
    }
    if (smallestEigenvalue(impedance.admittancePlus) <= 0.0 ||
        smallestEigenvalue(impedance.admittanceMinus) <= 0.0 ||
        smallestEigenvalue(impedance.eta) <= 0.0) {
      message << "fault impedance is not positive definite";
      return message.str();
    }
  }

  return std::nullopt;
}

} // namespace seissol::initializer::model

#endif // SEISSOL_SRC_INITIALIZER_MODEL_DYNAMICRUPTUREIMPEDANCE_H_
