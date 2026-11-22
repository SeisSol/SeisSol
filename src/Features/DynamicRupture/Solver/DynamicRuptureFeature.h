// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_DYNAMICRUPTUREFEATURE_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_DYNAMICRUPTUREFEATURE_H_

#include "Core/Kernels/GlobalDataStore.h"
#include "Core/Memory/Descriptor/DynamicRupture.h"
#include "Core/Run/Context.h"
#include "Core/TimeStepping/Compute/FaceClusterFactory.h"
#include "Core/TimeStepping/Feature.h"
#include "Features/DynamicRupture/Output/OutputManager.h"
#include "Features/DynamicRupture/Solver/DynamicRuptureClusterFactory.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolver.h"

namespace seissol::dr {

/**
 * The dynamic rupture as a feature of a run: the clusters of the faces of the fault, with the
 * friction law and the outputs on the fault. A run without a fault has no faces, so no clusters.
 */
class DynamicRuptureFeature : public solver::Feature {
  public:
  /// see DynamicRuptureClusterFactory
  DynamicRuptureFeature(RunContext& run,
                        DynamicRupture::Storage& storage,
                        GlobalDataStore& globalData,
                        const friction_law::FrictionSolverFactory& frictionLaw,
                        const friction_law::FrictionSolverFactory& frictionLawDevice,
                        output::OutputManager* outputManager,
                        double endTime)
      : clusters_(
            run, storage, globalData, frictionLaw, frictionLawDevice, outputManager, endTime) {}

  [[nodiscard]] const char* name() const override { return "dynamic rupture"; }

  [[nodiscard]] solver::FaceClusterFactory* faceClusters() override { return &clusters_; }

  private:
  DynamicRuptureClusterFactory clusters_;
};

} // namespace seissol::dr

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_DYNAMICRUPTUREFEATURE_H_
