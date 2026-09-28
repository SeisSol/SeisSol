// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_MODEL_DYNAMICRUPTUREIMPEDANCE_H_
#define SEISSOL_SRC_INITIALIZER_MODEL_DYNAMICRUPTUREIMPEDANCE_H_

#include "Equations/Impedance.h" // IWYU pragma: keep
#include "Equations/ImpedanceBase.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <optional>
#include <sstream>
#include <string>

namespace seissol::initializer::model {

/**
 * All impedance-derived quantities of a single fault face, in fault-local
 * coordinates (n, s, d [, fluid]). The admittances Y of the two sides come from
 * seissol::model::ImpedanceCompute (see there for why the stored quantity is an
 * admittance and not an impedance):
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
template <typename MaterialT>
struct FaultImpedance {
  using Matrix = typename seissol::model::ImpedanceCompute<MaterialT>::Matrix;
  using LateralMatrix = typename seissol::model::ImpedanceCompute<MaterialT>::LateralMatrix;

  Matrix admittancePlus;
  Matrix admittanceMinus;
  Matrix eta;
  Matrix bPlus;
  Matrix bMinus;
  /// only the plus side, matching the convention of ReceiverOutput::computeLocalStresses
  LateralMatrix lateralStressPlus;
};

template <typename MaterialT>
FaultImpedance<MaterialT>
    assembleFaultImpedance(const typename FaultImpedance<MaterialT>::Matrix& admittancePlus,
                           const typename FaultImpedance<MaterialT>::Matrix& admittanceMinus) {
  FaultImpedance<MaterialT> impedance;
  impedance.admittancePlus = admittancePlus;
  impedance.admittanceMinus = admittanceMinus;
  impedance.eta = (admittancePlus + admittanceMinus).inverse();
  impedance.bPlus = impedance.eta * admittancePlus;
  impedance.bMinus = impedance.eta * admittanceMinus;
  return impedance;
}

/// Both sides at once, materials given in the fault-local frame.
template <typename MaterialT>
FaultImpedance<MaterialT> computeFaultImpedance(const MaterialT& plusLocal,
                                                const MaterialT& minusLocal) {
  using LateralMatrix = typename FaultImpedance<MaterialT>::LateralMatrix;
  LateralMatrix lateralStressPlus = LateralMatrix::Zero();
  const auto admittancePlus = seissol::model::computeAdmittance(plusLocal, &lateralStressPlus);
  auto impedance = assembleFaultImpedance<MaterialT>(admittancePlus,
                                                     seissol::model::computeAdmittance(minusLocal));
  impedance.lateralStressPlus = lateralStressPlus;
  return impedance;
}

/**
 * Checks the invariants every fault impedance has to satisfy. Returns
 * std::nullopt on success, otherwise a human readable description of the first
 * violation.
 *
 * Self-adjointness is checked against the signature matrix of the material, so it also applies to
 * poroelastic, where the stored fourth traction component is +p while its energy conjugate partner
 * is -p. That part costs an eigendecomposition per call, everything before it a few matrix
 * products -- hence the switch, so that the cheap half can run in release builds as well.
 */
template <typename MaterialT>
std::optional<std::string> checkFaultImpedance(const FaultImpedance<MaterialT>& impedance,
                                               bool expectSelfAdjoint = true,
                                               double tolerance = 1e-9) {
  using Matrix = typename FaultImpedance<MaterialT>::Matrix;
  const Matrix signature = seissol::model::ImpedanceCompute<MaterialT>::signature();

  const auto relIdentityError = [](const Matrix& matrix) {
    return (matrix - Matrix::Identity()).cwiseAbs().maxCoeff();
  };
  const auto relSymmetryError = [&signature](const Matrix& matrix) {
    const Matrix adjusted = matrix * signature;
    const double scale = adjusted.cwiseAbs().maxCoeff();
    return (scale == 0.0) ? 0.0 : (adjusted - adjusted.transpose()).cwiseAbs().maxCoeff() / scale;
  };
  const auto smallestEigenvalue = [&signature](const Matrix& matrix) {
    const Matrix adjusted = matrix * signature;
    const Matrix symmetric = 0.5 * (adjusted + adjusted.transpose());
    const Eigen::SelfAdjointEigenSolver<Matrix> solver(symmetric, Eigen::EigenvaluesOnly);
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
      message << "fault impedance is not self-adjoint (w.r.t. the signature matrix), max relative "
                 "error "
              << symmetryError;
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
