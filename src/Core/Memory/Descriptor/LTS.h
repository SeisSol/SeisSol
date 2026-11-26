// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_CORE_MEMORY_DESCRIPTOR_LTS_H_
#define SEISSOL_SRC_CORE_MEMORY_DESCRIPTOR_LTS_H_

#include "Alignment.h"
#include "Core/Common/ConfigDispatch.h"
#include "Core/Common/Constants.h"
#include "Core/Memory/CellLocalInformation.h"
#include "Core/Memory/Tree/Backmap.h"
#include "Core/Memory/Tree/LTSTree.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Parallel/Helper.h"

/**
 * The storage of the cells, the layer-tree storage (LTS): its layers, by time cluster and halo
 * type, and the variables of the cells the core keeps itself. The wave propagation adds its
 * variables to the namespace and to the storage (see WavePropagation/Solver/Storage.h), and so do
 * the features (e.g. the plasticity, see solver::Feature::addStorage()).
 */
// (spelled as the storage of the cells always was, LTS, unlike other namespaces)
namespace seissol::LTS { // NOLINT(readability-identifier-naming)

enum class AllocationPreset {
  Global,
  Timedofs,
  Constant,
  Dofs,
  TimedofsConstant,
  ConstantShared,
  Timebucket,
  Plasticity,
  PlasticityData
};

/// The variables hold the data of every configuration; their sizes, which the placement
/// depends on, grow with the order.
inline auto allocationModeWP(AllocationPreset preset) {
  constexpr auto MaxOrder = maxBuiltConvergenceOrder();
  using namespace seissol::initializer;
  if constexpr (!isDeviceOn()) {
    switch (preset) {
    case AllocationPreset::Global:
      [[fallthrough]];
    case AllocationPreset::TimedofsConstant:
      return AllocationMode::HostOnlyHBM;
    case AllocationPreset::Plasticity:
      return AllocationMode::HostOnly;
    case AllocationPreset::PlasticityData:
      return AllocationMode::HostOnly;
    case AllocationPreset::Timebucket:
      [[fallthrough]];
    case AllocationPreset::Timedofs:
      return (MaxOrder <= 7 ? AllocationMode::HostOnlyHBM : AllocationMode::HostOnly);
    case AllocationPreset::Constant:
      [[fallthrough]];
    case AllocationPreset::ConstantShared:
      return (MaxOrder <= 4 ? AllocationMode::HostOnlyHBM : AllocationMode::HostOnly);
    case AllocationPreset::Dofs:
      return (MaxOrder <= 3 ? AllocationMode::HostOnlyHBM : AllocationMode::HostOnly);
    default:
      return AllocationMode::HostOnly;
    }
  } else {
    const auto modeMaybeCompress = useDeviceL2Compress() ? AllocationMode::HostDeviceCompress
                                                         : AllocationMode::HostDeviceSplit;
    const auto modeMaybeCompressPinned = useDeviceL2Compress()
                                             ? AllocationMode::HostDeviceCompressPinned
                                             : AllocationMode::HostDeviceSplitPinned;

    switch (preset) {
    case AllocationPreset::Global:
      [[fallthrough]];
    case AllocationPreset::Constant:
      [[fallthrough]];
    case AllocationPreset::TimedofsConstant:
      return AllocationMode::HostOnly;
    case AllocationPreset::Dofs:
      [[fallthrough]];
    case AllocationPreset::PlasticityData:
      return useUSM() ? AllocationMode::HostDeviceUnified : modeMaybeCompressPinned;
    case AllocationPreset::Timebucket:
      return useMPIUSM() ? AllocationMode::HostDeviceUnified : modeMaybeCompress;
    default:
      return useUSM() ? AllocationMode::HostDeviceUnified : modeMaybeCompress;
    }
  }
}

// the local information of a cell (its faces, its neighbors), and where it comes from
struct CellInformation : public initializer::Variable<CellLocalInformation> {};
struct SecondaryInformation : public initializer::Variable<SecondaryCellLocalInformation> {};

using initializer::Layer;
using initializer::Storage;
template <typename Cfg>
using Ref = initializer::CellRef<Cfg>;
using Backmap = initializer::StorageBackmap<Cell::NumFaces>;

/// Adds the variables of the cells the core keeps to `storage`, before it gets allocated; the wave
/// propagation and the features add theirs (see addWavePropagation()).
inline void addTo(Storage& storage) {
  using namespace initializer;
  storage.add<CellInformation>(
      LayerMask(), Alignment, allocationModeWP(AllocationPreset::Constant), true);
  storage.add<SecondaryInformation>(LayerMask(), Alignment, AllocationMode::HostOnly, true);
}

} // namespace seissol::LTS

#endif // SEISSOL_SRC_CORE_MEMORY_DESCRIPTOR_LTS_H_
