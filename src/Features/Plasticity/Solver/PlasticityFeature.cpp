// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "PlasticityFeature.h"

#include "Core/Common/ConfigDispatch.h"
#include "Core/Kernels/GlobalDataStore.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Run/Context.h"
#include "Core/TimeStepping/TimeManager.h"
#include "Features/Plasticity/PlasticityParameters.h"
#include "Features/Plasticity/Solver/PlasticityStage.h"
#include "Features/Plasticity/Solver/Storage.h"

#include <memory>
#include <utility>

namespace seissol::kernels {

PlasticityFeature::PlasticityFeature(RunContext& run,
                                     initializer::parameters::PlasticityParameters parameters,
                                     GlobalDataStore& globalData)
    : run_(run), parameters_(std::move(parameters)), globalData_(globalData) {}

void PlasticityFeature::addStorage(LTS::Storage& cells) {
  LTSPlasticity::addTo(cells, parameters_.enabled);
}

void PlasticityFeature::addStages(seissol::solver::TimeManager& timeManager, LTS::Storage& cells) {
  if (!parameters_.enabled) {
    return;
  }
  for (auto& layer : cells.leaves(Ghost)) {
    dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
      using Cfg = decltype(cfg);
      timeManager.addCellStage(layer.id(),
                               std::make_unique<PlasticityStage<Cfg>>(
                                   run_, parameters_.tv, globalData_.get<Cfg>(), layer));
    });
  }
}

} // namespace seissol::kernels
