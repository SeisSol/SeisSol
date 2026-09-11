// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_STP_SETUP_H_
#define SEISSOL_SRC_KERNELS_STP_SETUP_H_

#include "GeneratedCode/init.h"
#include "Kernels/STP/Solver.h"
#include "Model/Common.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cstddef>
#include <yateto.h>

namespace seissol::model {

/// True if the predictor has to factorise this row separately.
template <typename MaterialT>
constexpr bool isStiffRow(std::size_t quantity) {
  for (const auto& row : MaterialT::StiffSourceRows) {
    if (row.quantity == quantity) {
      return true;
    }
  }
  return false;
}

template <typename Tview>
inline void calcZinv(yateto::DenseTensorView<2, real, unsigned>& zInv,
                     const Tview& sourceMatrix,
                     size_t quantity,
                     bool isStiff,
                     double timeStepWidth) {
  using Matrix = Eigen::Matrix<real, ConvergenceOrder, ConvergenceOrder>;
  using Vector = Eigen::Matrix<real, ConvergenceOrder, 1>;

  Matrix matZ{init::Z::Values};
  // Only a stiff row carries a diagonal source entry. The check is not
  // cosmetic: for every other row the source matrix has no entry at
  // (quantity, quantity), so the lookup itself would be out of pattern.
  if (isStiff) {
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
template <typename MaterialT, size_t Istart, size_t Iend, typename Tview>
struct ZInvInitializer {
  ZInvInitializer(real* zInvData, const Tview& sourceMatrix, real timeStepWidth) {
    auto zInv = init::Zinv::view<Istart>::create(zInvData);
    calcZinv(zInv, sourceMatrix, Istart, isStiffRow<MaterialT>(Istart), timeStepWidth);
    if constexpr (Istart < Iend - 1) {
      auto* nextZInvData = zInvData + init::Zinv::size(Istart);
      ZInvInitializer<MaterialT, Istart + 1, Iend, Tview>(
          nextZInvData, sourceMatrix, timeStepWidth);
    }
  };
};

/**
 * The space-time predictor factorises the stiff rows of the source term
 * separately, so the per-cell data holds the inverses that go with them and
 * the off-diagonal entries they feed back. Which rows those are comes from the
 * material.
 */
template <typename MaterialT>
struct SolverSetup<kernels::solver::stp::Solver, MaterialT>
    : public SolverSetupDefaults<kernels::solver::stp::Solver, MaterialT> {
  static void initializeSpecificLocalData(const MaterialT& material,
                                          double timeStepWidth,
                                          typename MaterialT::Solver::LocalData* localData) {
    auto sourceMatrix = init::ET::view::create(localData->sourceMatrix);
    sourceMatrix.setZero();
    MaterialSetup<MaterialT>::getTransposedSourceCoefficientTensor(material, sourceMatrix);

    ZInvInitializer<MaterialT, 0, MaterialT::NumQuantities, decltype(sourceMatrix)>(
        localData->Zinv, sourceMatrix, timeStepWidth);

    std::fill(localData->G, localData->G + MaterialT::NumQuantities, 0.0);
    for (const auto& row : MaterialT::StiffSourceRows) {
      localData->G[row.quantity] = sourceMatrix(row.quantity, row.target);
    }

    localData->typicalTimeStepWidth = timeStepWidth;
  }
};

} // namespace seissol::model

#endif // SEISSOL_SRC_KERNELS_STP_SETUP_H_
