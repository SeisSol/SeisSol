// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_NONLINEARCK_SETUP_H_
#define SEISSOL_SRC_KERNELS_NONLINEARCK_SETUP_H_

// IWYU pragma: begin_exports

#include "Kernels/NonLinearCK/Local.h"
#include "Kernels/NonLinearCK/Neighbor.h"
#include "Kernels/NonLinearCK/Time.h"

// IWYU pragma: end_exports

#include "GeneratedCode/init.h"
#include "GeneratedCode/quantities.h"
#include "Kernels/NonLinearCK/Solver.h"
#include "Model/Common.h"

#include <cstddef>

namespace seissol::model {

/**
 * The initial strain is a material parameter that the kernels read as a
 * tensor, so it is converted once, here, rather than at every timestep: the
 * material holds it in double and the kernels work in the solver's precision.
 */
template <typename MaterialT>
struct SolverSetup<kernels::solver::nonlinearck::Solver, MaterialT>
    : public SolverSetupDefaults<kernels::solver::nonlinearck::Solver, MaterialT> {
  // The transported quantities are laid out on both sides independently, and
  // a ghost rule rotates them, so a disagreement would be a wrong reflection
  // rather than a compile error. This makes it one.
  static_assert(detail::kindsMatch(MaterialT::TransportGroups,
                                   generated::TransportRotationGroupKinds),
                "The transport layout of the material and of the generated code disagree.");

  static void initializeSpecificLocalData(const MaterialT& material,
                                          double /*timeStepWidth*/,
                                          typename MaterialT::Solver::LocalData* localData) {
    auto epsInit = init::epsInit::view::create(localData->epsInit);
    epsInit.setZero();
    epsInit(0) = material.epsInitXX;
    epsInit(1) = material.epsInitYY;
    epsInit(2) = material.epsInitZZ;
    epsInit(3) = material.epsInitXY;
    epsInit(4) = material.epsInitYZ;
    epsInit(5) = material.epsInitXZ;

    localData->maxWaveSpeedBound = static_cast<real>(material.getMaxWaveSpeed());

    auto& parameters = localData->parameters;
    parameters.rhoInv = static_cast<real>(1.0 / material.rho);
    parameters.lambda0 = static_cast<real>(material.lambda0);
    parameters.mu0 = static_cast<real>(material.mu0);
    parameters.gammaR = static_cast<real>(material.gammaR);
    parameters.xi0 = static_cast<real>(material.xi0);
    parameters.damageRate = static_cast<real>(material.damageRate);
    parameters.breakageRate = static_cast<real>(material.breakageRate);
    parameters.healingRate = static_cast<real>(material.healingRate);
    parameters.betaAlpha = static_cast<real>(material.betaAlpha);
    for (std::size_t i = 0; i < parameters.aB.size(); ++i) {
      parameters.aB[i] = static_cast<real>(material.aB.at(i));
    }
  }

  static void
      initializeSpecificNeighborData(const MaterialT& material,
                                     typename MaterialT::Solver::NeighborData* neighborData) {
    neighborData->maxWaveSpeedBound.fill(static_cast<real>(material.getMaxWaveSpeed()));
  }
};

} // namespace seissol::model

#endif // SEISSOL_SRC_KERNELS_NONLINEARCK_SETUP_H_
