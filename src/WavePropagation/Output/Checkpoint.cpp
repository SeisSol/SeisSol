// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Checkpoint.h"

#include "Core/Common/ConfigDispatch.h"
#include "Core/Common/ConfigRegistry.h"
#include "Core/IO/Instance/Checkpoint/CheckpointManager.h"
#include "Core/Kernels/Common.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Tree/Layer.h"
#include "GeneratedCode/tensor.h"
#include "WavePropagation/Output/SurfaceStorage.h"
#include "WavePropagation/Solver/Storage.h"

namespace seissol::solver {

void registerCellCheckpoint(io::instance::checkpoint::CheckpointManager& manager,
                            LTS::Storage& storage,
                            ConfigId config) {
  manager.registerData<LTS::Dofs>("dofs", storage);
  dispatchConfig(config, [&](auto cfg) {
    if constexpr (kernels::size<tensor::Qane<decltype(cfg)>>() > 0) {
      manager.registerData<LTS::DofsAne>("dofsAne", storage);
    }
  });
  // the time integrals of the unknowns, if the output integrates them
  if (storage.info<LTS::Integrals>().mask == initializer::LayerMask(Ghost)) {
    manager.registerData<LTS::Integrals>("integrals", storage);
  }
}

void registerSurfaceCheckpoint(io::instance::checkpoint::CheckpointManager& manager,
                               SurfaceLTS::Storage& storage) {
  manager.registerData<SurfaceLTS::DisplacementDofs>("displacementDofs", storage);
}

} // namespace seissol::solver
