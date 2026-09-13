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
#include <string>
#include <string_view>
#include <utils/logger.h>

namespace seissol::model {

/**
 * The initial strain is a material parameter that the kernels read as a
 * tensor, so it is converted once, here, rather than at every timestep: the
 * material holds it in double and the kernels work in the solver's precision.
 */
template <typename MaterialT>
struct SolverSetup<kernels::solver::nonlinearck::Solver, MaterialT>
    : public SolverSetupDefaults<kernels::solver::nonlinearck::Solver, MaterialT> {
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

    // Filled by name against the order the codegen chose, so the two cannot
    // drift apart: a parameter in the wrong slot would be a different
    // material, silently.
    const auto set = [&](std::string_view name, double value) {
      for (std::size_t i = 0; i < generated::MaterialParameterNames.size(); ++i) {
        if (generated::MaterialParameterNames[i] == name) {
          localData->parameters[i] = static_cast<real>(value);
          return;
        }
      }
      logError() << "The kernels do not read a material parameter called" << name.data();
    };

    set("rhoInv", 1.0 / material.rho);
    set("lambda0", material.lambda0);
    set("mu0", material.mu0);
    set("gammaR", material.gammaR);
    set("xi0", material.xi0);
    set("damageRate", material.damageRate);
    set("breakageRate", material.breakageRate);
    set("healingRate", material.healingRate);
    set("betaAlpha", material.betaAlpha);
    for (std::size_t i = 0; i < material.aB.size(); ++i) {
      set("aB" + std::to_string(i), material.aB.at(i));
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
