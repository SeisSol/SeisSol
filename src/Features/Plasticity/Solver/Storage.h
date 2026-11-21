// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_PLASTICITY_SOLVER_STORAGE_H_
#define SEISSOL_SRC_FEATURES_PLASTICITY_SOLVER_STORAGE_H_

#include "Alignment.h"
#include "Core/Common/Real.h"
#include "Core/IO/Instance/Checkpoint/CheckpointManager.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Parallel/Helper.h"
#include "Features/Plasticity/Physics/Plasticity.h"
#include "GeneratedCode/tensor.h"

namespace seissol {

/**
 * The data of the plasticity in the storage of the cells (see LTS), which the plasticity adds there
 * itself (see addTo()).
 */
struct LTSPlasticity {
  // the plastic strain of a cell and its accumulated plastic strain (eta), nodal
  template <typename Cfg>
  using PStrainArray =
      Real<Cfg>[tensor::QStressNodal<Cfg>::size() + tensor::QEtaNodal<Cfg>::size()];

  struct Plasticity : public initializer::VariantVariable<model::PlasticityData> {};
  struct PStrain : public initializer::VariantVariable<PStrainArray> {};

  // (on the device) whether the plastic cells of a layer yield, and their nodal stresses
  struct FlagScratch : public initializer::Scratchpad<unsigned> {};
  struct QStressNodalScratch : public initializer::VariantScratchpad<Real> {};

  /// Adds the data to the storage of the cells; without plasticity, for no layer.
  static void addTo(LTS::Storage& storage, bool plasticity) {
    using namespace initializer;
    const auto mask =
        plasticity ? LayerMask(Ghost) : LayerMask(Ghost) | LayerMask(Copy) | LayerMask(Interior);
    storage.add<Plasticity>(
        mask, Alignment, LTS::allocationModeWP(LTS::AllocationPreset::Plasticity), true);
    storage.add<PStrain>(
        mask, PagesizeHeap, LTS::allocationModeWP(LTS::AllocationPreset::PlasticityData));

    if constexpr (isDeviceOn()) {
      storage.add<FlagScratch>(LayerMask(), Alignment, AllocationMode::DeviceOnly);
      storage.add<QStressNodalScratch>(LayerMask(), Alignment, AllocationMode::DeviceOnly);
    }
  }

  /// The data to checkpoint, if there is plasticity.
  static void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                          LTS::Storage& storage) {
    // check plasticity usage over the layer mask (for now)
    if (storage.info<Plasticity>().mask == initializer::LayerMask(Ghost)) {
      manager.registerData<PStrain>("pstrain", storage);
    }
  }
};

} // namespace seissol

#endif // SEISSOL_SRC_FEATURES_PLASTICITY_SOLVER_STORAGE_H_
