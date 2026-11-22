// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "IntegrationFeature.h"

#include "Core/Common/ConfigDispatch.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/TimeStepping/TimeManager.h"
#include "WavePropagation/Output/IntegrationStage.h"

#include <memory>

namespace seissol::kernels {

void IntegrationFeature::addStages(seissol::solver::TimeManager& timeManager, LTS::Storage& cells) {
  if (!integrate_) {
    return;
  }
  for (auto& layer : cells.leaves(Ghost)) {
    dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
      using Cfg = decltype(cfg);
      timeManager.addCellStage(layer.id(), std::make_unique<IntegrationStage<Cfg>>(layer));
    });
  }
}

} // namespace seissol::kernels
