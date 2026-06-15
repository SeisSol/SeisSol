// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MEMORY_DESCRIPTOR_BOUNDARY_H_
#define SEISSOL_SRC_MEMORY_DESCRIPTOR_BOUNDARY_H_

#include "Alignment.h"
#include "IO/Instance/Checkpoint/CheckpointManager.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"
#include "Parallel/Helper.h"

namespace seissol {

inline auto allocationModeBoundary() {
  using namespace initializer;
  if constexpr (!isDeviceOn()) {
    return AllocationMode::HostOnly;
  } else {
    const auto modeMaybeCompress = useDeviceL2Compress() ? AllocationMode::HostDeviceCompress
                                                         : AllocationMode::HostDeviceSplit;
    return useUSM() ? AllocationMode::HostDeviceUnified : modeMaybeCompress;
  }
}

struct Boundary {
  struct FaceInformation : public initializer::Variable<BoundaryFaceInformation> {};

  struct BoundaryVarmap : public initializer::SpecificVarmap<FaceInformation> {};

  using Storage = initializer::Storage<BoundaryVarmap>;
  using Layer = initializer::Layer<BoundaryVarmap>;
  using Ref = initializer::Layer<BoundaryVarmap>::CellRef;
  using Backmap = initializer::StorageBackmap<1>;

  static void addTo(Storage& storage) {
    const auto mask = initializer::LayerMask(Ghost);
    storage.add<FaceInformation>(mask, Alignment, allocationModeBoundary());
  }
};

} // namespace seissol

#endif // SEISSOL_SRC_MEMORY_DESCRIPTOR_BOUNDARY_H_
