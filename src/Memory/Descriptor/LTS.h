// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_MEMORY_DESCRIPTOR_LTS_H_
#define SEISSOL_SRC_MEMORY_DESCRIPTOR_LTS_H_

#include "GeneratedCode/tensor.h"
#include "IO/Instance/Checkpoint/CheckpointManager.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"
#include "Model/Plasticity.h"
#include <Equations/Datastructures.h>
#include <Initializer/CellLocalInformation.h>
#include <Memory/Tree/Backmap.h>
#include <Parallel/Helper.h>

#ifdef ACL_DEVICE
#include "Parallel/Helper.h"
#endif

namespace seissol {

struct LTS {
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

  static auto allocationModeWP(AllocationPreset preset, int convergenceOrder = 6) {
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
        return (convergenceOrder <= 7 ? AllocationMode::HostOnlyHBM : AllocationMode::HostOnly);
      case AllocationPreset::Constant:
        [[fallthrough]];
      case AllocationPreset::ConstantShared:
        return (convergenceOrder <= 4 ? AllocationMode::HostOnlyHBM : AllocationMode::HostOnly);
      case AllocationPreset::Dofs:
        return (convergenceOrder <= 3 ? AllocationMode::HostOnlyHBM : AllocationMode::HostOnly);
      default:
        return AllocationMode::HostOnly;
      }
    } else {
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
        return useUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplitPinned;
      case AllocationPreset::Timebucket:
        return useMPIUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplit;
      default:
        return useUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplit;
      }
    }
  }

  // to prevent a clang-format bug
  template <typename Cfg>
  using RealPtr = Real<Cfg>*;

  struct Dofs : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = Real<Cfg>[tensor::Q<Cfg>::size()];
  };
  // size is zero if Qane is not defined
  struct DofsAne : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = Real<Cfg>[zeroLengthArrayHandler(kernels::size<tensor::Qane<Cfg>>())];
  };
  struct Buffers : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = RealPtr<Cfg>;
  };
  struct Derivatives : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = RealPtr<Cfg>;
  };
  struct CellInformation : public initializer::Variable<CellLocalInformation> {};
  struct SecondaryInformation : public initializer::Variable<SecondaryCellLocalInformation> {};
  struct FaceNeighbors : public initializer::Variable<void* [Cell::NumFaces]> {};
  struct LocalIntegration : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = LocalIntegrationData<Cfg>;
  };
  struct NeighboringIntegration : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = NeighboringIntegrationData<Cfg>;
  };
  struct MaterialData : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = model::MaterialTT<Cfg>;
  };
  struct Material : public initializer::Variable<CellMaterialData> {};
  struct Plasticity : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = seissol::model::PlasticityData<Cfg>;
  };
  struct DRMapping : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = CellDRMapping<Cfg>[Cell::NumFaces];
  };
  struct BoundaryMapping : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = CellBoundaryMapping<Cfg>[Cell::NumFaces];
  };
  struct PStrain : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = Real<Cfg>[tensor::QStress<Cfg>::size() + tensor::QEtaModal<Cfg>::size()];
  };
  struct FaceDisplacements : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = RealPtr<Cfg>[Cell::NumFaces];
  };
  struct BuffersDerivatives : public initializer::Bucket<void> {
    template <typename Cfg>
    using VariantType = Real<Cfg>;
  };

  struct BuffersDevice : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = RealPtr<Cfg>;
  };
  struct DerivativesDevice : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = RealPtr<Cfg>;
  };
  struct FaceNeighborsDevice : public initializer::Variable<void* [Cell::NumFaces]> {};
  struct FaceDisplacementsDevice : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = RealPtr<Cfg>[Cell::NumFaces];
  };
  struct DRMappingDevice : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = CellDRMapping<Cfg>[Cell::NumFaces];
  };
  struct BoundaryMappingDevice : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = CellBoundaryMapping<Cfg>[Cell::NumFaces];
  };

  struct ScratchBase : public initializer::Scratchpad<void> {
    template <typename Cfg>
    using VariantType = Real<Cfg>;
  };

  struct IntegratedDofsScratch : public ScratchBase {};
  struct DerivativesScratch : public ScratchBase {};
  struct NodalAvgDisplacements : public ScratchBase {};
  struct AnalyticScratch : public ScratchBase {};
  struct DerivativesExtScratch : public ScratchBase {};
  struct DerivativesAneScratch : public ScratchBase {};
  struct IDofsAneScratch : public ScratchBase {};
  struct DofsExtScratch : public ScratchBase {};

  struct FlagScratch : public initializer::Scratchpad<unsigned> {};
  struct PrevDofsScratch : public ScratchBase {};
  struct QEtaNodalScratch : public ScratchBase {};
  struct QStressNodalScratch : public ScratchBase {};

  struct RotateDisplacementToFaceNormalScratch : public ScratchBase {};
  struct RotateDisplacementToGlobalScratch : public ScratchBase {};
  struct RotatedFaceDisplacementScratch : public ScratchBase {};
  struct DofsFaceNodalScratch : public ScratchBase {};
  struct PrevCoefficientsScratch : public ScratchBase {};
  struct DofsFaceBoundaryNodalScratch : public ScratchBase {};

  struct Integrals : public initializer::VariantVariable<Real> {};

  struct LTSVarmap : public initializer::SpecificVarmap<Dofs,
                                                        DofsAne,
                                                        Buffers,
                                                        Derivatives,
                                                        CellInformation,
                                                        SecondaryInformation,
                                                        FaceNeighbors,
                                                        LocalIntegration,
                                                        NeighboringIntegration,
                                                        Material,
                                                        MaterialData,
                                                        Plasticity,
                                                        DRMapping,
                                                        BoundaryMapping,
                                                        PStrain,
                                                        FaceDisplacements,
                                                        BuffersDerivatives,
                                                        BuffersDevice,
                                                        DerivativesDevice,
                                                        FaceNeighborsDevice,
                                                        FaceDisplacementsDevice,
                                                        DRMappingDevice,
                                                        BoundaryMappingDevice,
                                                        IntegratedDofsScratch,
                                                        DerivativesScratch,
                                                        NodalAvgDisplacements,
                                                        AnalyticScratch,
                                                        DerivativesExtScratch,
                                                        DerivativesAneScratch,
                                                        IDofsAneScratch,
                                                        DofsExtScratch,
                                                        FlagScratch,
                                                        PrevDofsScratch,
                                                        QEtaNodalScratch,
                                                        QStressNodalScratch,
                                                        RotateDisplacementToFaceNormalScratch,
                                                        RotateDisplacementToGlobalScratch,
                                                        RotatedFaceDisplacementScratch,
                                                        DofsFaceNodalScratch,
                                                        PrevCoefficientsScratch,
                                                        DofsFaceBoundaryNodalScratch,
                                                        Integrals> {};

  using Storage = initializer::Storage<LTSVarmap>;
  using Layer = initializer::Layer<LTSVarmap>;
  template <typename Config>
  using Ref = initializer::Layer<LTSVarmap>::CellRef<Config>;
  using Backmap = initializer::StorageBackmap<Cell::NumFaces>;

  static void addTo(Storage& storage, bool usePlasticity) {
    using namespace initializer;
    LayerMask plasticityMask;
    if (usePlasticity) {
      plasticityMask = LayerMask(Ghost);
    } else {
      plasticityMask = LayerMask(Ghost) | LayerMask(Copy) | LayerMask(Interior);
    }

    storage.add<Dofs>(LayerMask(Ghost), PagesizeHeap, allocationModeWP(AllocationPreset::Dofs));

    storage.add<DofsAne>(LayerMask(Ghost), PagesizeHeap, allocationModeWP(AllocationPreset::Dofs));

    storage.add<Buffers>(
        LayerMask(), 1, allocationModeWP(AllocationPreset::TimedofsConstant), true);
    storage.add<Derivatives>(
        LayerMask(), 1, allocationModeWP(AllocationPreset::TimedofsConstant), true);
    storage.add<CellInformation>(
        LayerMask(), 1, allocationModeWP(AllocationPreset::Constant), true);
    storage.add<SecondaryInformation>(LayerMask(), 1, AllocationMode::HostOnly, true);
    storage.add<FaceNeighbors>(
        LayerMask(Ghost), 1, allocationModeWP(AllocationPreset::TimedofsConstant), true);
    storage.add<LocalIntegration>(
        LayerMask(Ghost), 1, allocationModeWP(AllocationPreset::ConstantShared), true);
    storage.add<NeighboringIntegration>(
        LayerMask(Ghost), 1, allocationModeWP(AllocationPreset::ConstantShared), true);
    storage.add<MaterialData>(LayerMask(), 1, AllocationMode::HostOnly, true);
    storage.add<Material>(LayerMask(Ghost), 1, AllocationMode::HostOnly, true);
    storage.add<Plasticity>(
        plasticityMask, 1, allocationModeWP(AllocationPreset::Plasticity), true);
    storage.add<DRMapping>(LayerMask(Ghost), 1, allocationModeWP(AllocationPreset::Constant), true);
    storage.add<BoundaryMapping>(
        LayerMask(Ghost), 1, allocationModeWP(AllocationPreset::Constant), true);
    storage.add<PStrain>(
        plasticityMask, PagesizeHeap, allocationModeWP(AllocationPreset::PlasticityData));
    storage.add<FaceDisplacements>(LayerMask(Ghost), PagesizeHeap, AllocationMode::HostOnly, true);

    // TODO(David): remove/rename "constant" flag (the data is temporary; and copying it for IO is
    // handled differently)
    storage.add<BuffersDerivatives>(
        LayerMask(), PagesizeHeap, allocationModeWP(AllocationPreset::Timebucket), true);

    storage.add<BuffersDevice>(LayerMask(), 1, AllocationMode::HostOnly, true);
    storage.add<DerivativesDevice>(LayerMask(), 1, AllocationMode::HostOnly, true);
    storage.add<FaceDisplacementsDevice>(LayerMask(Ghost), 1, AllocationMode::HostOnly, true);
    storage.add<FaceNeighborsDevice>(LayerMask(Ghost), 1, AllocationMode::HostOnly, true);
    storage.add<DRMappingDevice>(LayerMask(Ghost), 1, AllocationMode::HostOnly, true);
    storage.add<BoundaryMappingDevice>(LayerMask(Ghost), 1, AllocationMode::HostOnly, true);

    if constexpr (isDeviceOn()) {
      storage.add<DerivativesExtScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      storage.add<DerivativesAneScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      storage.add<IDofsAneScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      storage.add<DofsExtScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      storage.add<IntegratedDofsScratch>(LayerMask(), 1, AllocationMode::HostDeviceSplit);
      storage.add<DerivativesScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      storage.add<NodalAvgDisplacements>(LayerMask(), 1, AllocationMode::DeviceOnly);
      storage.add<AnalyticScratch>(LayerMask(), 1, AllocationMode::HostDevicePinned);

      storage.add<FlagScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      storage.add<PrevDofsScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      storage.add<QEtaNodalScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      storage.add<QStressNodalScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);

      storage.add<RotateDisplacementToFaceNormalScratch>(
          LayerMask(), 1, AllocationMode::DeviceOnly);
      storage.add<RotateDisplacementToGlobalScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      storage.add<RotatedFaceDisplacementScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      storage.add<DofsFaceNodalScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      storage.add<PrevCoefficientsScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      storage.add<DofsFaceBoundaryNodalScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
    }
  }

  static void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                          Storage& storage) {
    manager.registerData<Dofs>("dofs", storage);
    manager.registerData<DofsAne>("dofsAne", storage);
    // check plasticity usage over the layer mask (for now)
    if (storage.info<Plasticity>().mask == initializer::LayerMask(Ghost)) {
      manager.registerData<Plasticity>("pstrain", storage);
    }
  }
};

} // namespace seissol

#endif // SEISSOL_SRC_MEMORY_DESCRIPTOR_LTS_H_
