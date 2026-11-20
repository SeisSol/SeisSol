// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_MEMORY_DESCRIPTOR_BOUNDARY_H_
#define SEISSOL_SRC_CORE_MEMORY_DESCRIPTOR_BOUNDARY_H_

#include "Alignment.h"
#include "Core/IO/Instance/Checkpoint/CheckpointManager.h"
#include "Core/Memory/Tree/LTSTree.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Parallel/Helper.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"

namespace seissol {

inline auto allocationModeBoundary() {
  using namespace initializer;
  if constexpr (!isDeviceOn()) {
    return AllocationMode::HostOnly;
  } else {
    // pinned on the host, so that the copies of the data for the outputs do not block (see
    // HostSnapshot)
    const auto modeMaybeCompress = useDeviceL2Compress() ? AllocationMode::HostDeviceCompressPinned
                                                         : AllocationMode::HostDeviceSplitPinned;
    return useUSM() ? AllocationMode::HostDeviceUnified : modeMaybeCompress;
  }
}

struct Boundary {
  struct FaceInformation : public initializer::VariantVariable<BoundaryFaceInformation> {};

  using Storage = initializer::Storage;
  using Layer = initializer::Layer;
  template <typename Cfg>
  using Ref = initializer::CellRef<Cfg>;
  using Backmap = initializer::StorageBackmap<1>;

  static void addTo(Storage& storage) {
    const auto mask = initializer::LayerMask(Ghost);
    storage.add<FaceInformation>(mask, Alignment, allocationModeBoundary());
  }
};

} // namespace seissol

#endif // SEISSOL_SRC_CORE_MEMORY_DESCRIPTOR_BOUNDARY_H_
