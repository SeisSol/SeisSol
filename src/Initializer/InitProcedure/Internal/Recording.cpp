// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Recording.h"

#include "Initializer/BatchRecorders/Recorders.h"
#include "Kernels/Common.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"

namespace seissol::initializer::internal {

void setupRecorders(LTS::Storage& ltsStorage,
                    DynamicRupture::Storage& drStorage,
                    bool usePlasticity) {
  // only run for GPUs
  if constexpr (isDeviceOn()) {
    recording::CompositeRecorder<LTS::LTSVarmap> recorder;
    recorder.addRecorder(new recording::LocalIntegrationRecorder);
    recorder.addRecorder(new recording::NeighIntegrationRecorder);

    if (usePlasticity) {
      recorder.addRecorder(new recording::PlasticityRecorder);
    }

    for (auto& layer : ltsStorage.leaves(Ghost)) {
      recorder.record(layer);
    }

    recording::CompositeRecorder<DynamicRupture::DynrupVarmap> drRecorder;
    drRecorder.addRecorder(new recording::DynamicRuptureRecorder);
    for (auto& layer : drStorage.leaves(Ghost)) {
      drRecorder.record(layer);
    }
  }
}

} // namespace seissol::initializer::internal
