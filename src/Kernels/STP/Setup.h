// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_STP_SETUP_H_
#define SEISSOL_SRC_KERNELS_STP_SETUP_H_

#include "Equations/poroelastic/Model/Helper.h"
#include "GeneratedCode/init.h"
#include "Kernels/STP/Solver.h"
#include "Model/Common.h"

#include <cstddef>

namespace seissol::model {

/**
 * The space-time predictor peels the stiff rows of the source term into a
 * factorisation of its own, so the per-cell data holds those rows and the
 * inverses that go with them. Which rows they are is still written out here;
 * it ought to come from the material.
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

    ZInvInitializer<0, PoroElasticMaterial::NumQuantities, decltype(sourceMatrix)>(
        localData->Zinv, sourceMatrix, timeStepWidth);
    std::fill(localData->G, localData->G + PoroElasticMaterial::NumQuantities, 0.0);
    localData->G[10] = sourceMatrix(10, 6);
    localData->G[11] = sourceMatrix(11, 7);
    localData->G[12] = sourceMatrix(12, 8);

    localData->typicalTimeStepWidth = timeStepWidth;
  }
};

} // namespace seissol::model

#endif // SEISSOL_SRC_KERNELS_STP_SETUP_H_
