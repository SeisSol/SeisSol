// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_PLASTICITY_SOLVER_PLASTICITYFEATURE_H_
#define SEISSOL_SRC_FEATURES_PLASTICITY_SOLVER_PLASTICITYFEATURE_H_

#include "Core/Kernels/GlobalDataStore.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Run/Context.h"
#include "Core/TimeStepping/Feature.h"
#include "Features/Plasticity/PlasticityParameters.h"

namespace seissol::kernels {

/**
 * The plasticity as a feature of a run: its variables in the storage of the cells (without memory
 * if the run has no plasticity), and a cell stage per cluster that lets the cells yield after
 * each correction.
 */
class PlasticityFeature : public seissol::solver::Feature {
  public:
  PlasticityFeature(RunContext& run,
                    initializer::parameters::PlasticityParameters parameters,
                    GlobalDataStore& globalData);

  [[nodiscard]] const char* name() const override { return "plasticity"; }

  void addStorage(LTS::Storage& cells) override;

  void addStages(seissol::solver::TimeManager& timeManager, LTS::Storage& cells) override;

  private:
  RunContext& run_;
  initializer::parameters::PlasticityParameters parameters_;
  GlobalDataStore& globalData_;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_FEATURES_PLASTICITY_SOLVER_PLASTICITYFEATURE_H_
