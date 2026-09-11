// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_LINEARCK_SETUP_H_
#define SEISSOL_SRC_KERNELS_LINEARCK_SETUP_H_

#include "GeneratedCode/init.h"
#include "Kernels/LinearCK/Solver.h"
#include "Model/Common.h"

namespace seissol::model {

/**
 * The memory variables, where there are any, sit in Q alongside the elastic
 * quantities. The whole system is therefore one operator, and the source term
 * is one matrix over all of it -- which is why the plane wave operator needs
 * nothing beyond the defaults.
 */
template <typename MaterialT>
struct SolverSetup<kernels::solver::linearck::Solver, MaterialT>
    : public SolverSetupDefaults<kernels::solver::linearck::Solver, MaterialT> {
  static void initializeSpecificLocalData(const MaterialT& material,
                                          double /*timeStepWidth*/,
                                          typename MaterialT::Solver::LocalData* localData) {
    auto sourceMatrix = init::ET::view::create(localData->sourceMatrix);
    sourceMatrix.setZero();
    MaterialSetup<MaterialT>::getTransposedSourceCoefficientTensor(material, sourceMatrix);
  }
};

} // namespace seissol::model

#endif // SEISSOL_SRC_KERNELS_LINEARCK_SETUP_H_
