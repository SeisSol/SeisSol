// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MEMORY_DESCRIPTOR_BOUNDARY_H_
#define SEISSOL_SRC_MEMORY_DESCRIPTOR_BOUNDARY_H_

#include "IO/Instance/Checkpoint/CheckpointManager.h"
#include "Initializer/Typedefs.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"
#include "Parallel/Helper.h"
#include <Kernels/Common.h>

namespace seissol {

inline auto allocationModeBoundary() {
  using namespace initializer;
  if constexpr (!isDeviceOn()) {
    return AllocationMode::HostOnly;
  } else {
    return useUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplit;
  }
}

struct Boundary {
  struct FaceInformation : public initializer::Variable<BoundaryFaceInformation> {};

  void addTo(initializer::LTSTree& tree) {
    const auto mask = initializer::LayerMask(Ghost);
    tree.add<FaceInformation>(mask, 1, allocationModeBoundary());
  }
};

} // namespace seissol

#endif // SEISSOL_SRC_MEMORY_DESCRIPTOR_BOUNDARY_H_
