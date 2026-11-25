// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Recording.h"

#include "Common/ConfigDispatch.h"
#include "Initializer/BatchRecorders/Recorders.h"
#include "Kernels/Common.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"

namespace seissol::initializer::internal {

void setupRecorders(LTS::Storage& ltsStorage,
                    DynamicRupture::Storage& drStorage,
                    bool usePlasticity,
                    double g) {
  // only run for GPUs
  if constexpr (isDeviceOn()) {
    // every layer is recorded in its configuration
    for (auto& layer : ltsStorage.leaves(Ghost)) {
      dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
        using Cfg = decltype(cfg);
        recording::CompositeRecorder<LTS::LTSVarmap> recorder;
        recorder.addRecorder(new recording::LocalIntegrationRecorder<Cfg>(g));
        recorder.addRecorder(new recording::NeighIntegrationRecorder<Cfg>());
        if (usePlasticity) {
          recorder.addRecorder(new recording::PlasticityRecorder<Cfg>());
        }
        recorder.record(layer);
      });
    }

    for (auto& layer : drStorage.leaves(Ghost)) {
      dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
        using Cfg = decltype(cfg);
        recording::CompositeRecorder<DynamicRupture::DynrupVarmap> drRecorder;
        drRecorder.addRecorder(new recording::DynamicRuptureRecorder<Cfg>());
        drRecorder.record(layer);
      });
    }
  }
}

} // namespace seissol::initializer::internal
