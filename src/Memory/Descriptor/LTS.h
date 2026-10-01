// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_MEMORY_DESCRIPTOR_LTS_H_
#define SEISSOL_SRC_MEMORY_DESCRIPTOR_LTS_H_

#include "Alignment.h"
#include "Common/Real.h"
#include "Config.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/tensor.h"
#include "IO/Instance/Checkpoint/CheckpointManager.h"
#include "Initializer/CellLocalInformation.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Memory/Tree/Backmap.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"
#include "Model/Plasticity.h"
#include "Parallel/Helper.h"
#include "Solver/Settings.h"

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

  static auto allocationModeWP(AllocationPreset preset,
                               int convergenceOrder = seissol::ConvergenceOrder) {
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

  // The unknowns of a cell and what is derived from them are held in the reals and the layout of
  // the configuration of their layer.
  template <typename Cfg>
  using DofsArray = Real<Cfg>[tensor::Q<Cfg>::size()];
  // empty if the configuration has no Qane
  template <typename Cfg>
  using DofsAneArray = Real<Cfg>[zeroGuard(kernels::size<tensor::Qane<Cfg>>())];
  template <typename Cfg>
  using PStrainArray =
      Real<Cfg>[tensor::QStressNodal<Cfg>::size() + tensor::QEtaNodal<Cfg>::size()];
  template <typename Cfg>
  using RealPtr = Real<Cfg>*;
  template <typename Cfg>
  using FaceRealPtrs = std::array<Real<Cfg>*, Cell::NumFaces>;

  struct Dofs : public initializer::VariantVariable<DofsArray> {};
  struct DofsHalo : public initializer::VariantVariable<DofsArray> {};
  struct DofsAne : public initializer::VariantVariable<DofsAneArray> {};
  struct StepIntegrals : public initializer::VariantVariable<RealPtr> {};
  struct AccumulatedIntegrals : public initializer::VariantVariable<RealPtr> {};
  struct Derivatives : public initializer::VariantVariable<RealPtr> {};
  struct CellInformation : public initializer::Variable<CellLocalInformation> {};
  struct SecondaryInformation : public initializer::Variable<SecondaryCellLocalInformation> {};
  struct FaceNeighbors : public initializer::Variable<std::array<real*, Cell::NumFaces>> {};
  struct LocalIntegration : public initializer::Variable<LocalIntegrationData> {};
  struct NeighboringIntegration : public initializer::Variable<NeighboringIntegrationData> {};
  struct MaterialData : public initializer::Variable<model::MaterialT> {};
  struct Material : public initializer::Variable<CellMaterialData> {};
  struct Plasticity : public initializer::Variable<seissol::model::PlasticityData> {};
  struct DRMapping : public initializer::Variable<std::array<CellDRMapping, Cell::NumFaces>> {};
  struct BoundaryMapping
      : public initializer::Variable<std::array<CellBoundaryMapping, Cell::NumFaces>> {};
  struct PStrain : public initializer::VariantVariable<PStrainArray> {};
  struct FaceDisplacements : public initializer::VariantVariable<FaceRealPtrs> {};
  struct Buffers : public initializer::VariantBucket<Real> {};

  struct StepIntegralsDevice : public initializer::VariantVariable<RealPtr> {};
  struct AccumulatedIntegralsDevice : public initializer::VariantVariable<RealPtr> {};
  struct DerivativesDevice : public initializer::VariantVariable<RealPtr> {};
  struct FaceNeighborsDevice : public initializer::Variable<std::array<real*, Cell::NumFaces>> {};
  struct FaceDisplacementsDevice : public initializer::VariantVariable<FaceRealPtrs> {};
  struct DRMappingDevice : public initializer::Variable<std::array<CellDRMapping, Cell::NumFaces>> {
  };
  struct BoundaryMappingDevice
      : public initializer::Variable<std::array<CellBoundaryMapping, Cell::NumFaces>> {};

  struct EnergyData : public initializer::Variable<model::MaterialT::EnergyData> {};

  struct IntegratedDofsScratch : public initializer::VariantScratchpad<Real> {};
  struct DerivativesScratch : public initializer::VariantScratchpad<Real> {};
  struct NodalAvgDisplacements : public initializer::VariantScratchpad<Real> {};
  struct AnalyticScratch : public initializer::VariantScratchpad<Real> {};
  struct DerivativesExtScratch : public initializer::VariantScratchpad<Real> {};
  struct DerivativesAneScratch : public initializer::VariantScratchpad<Real> {};
  struct IDofsAneScratch : public initializer::VariantScratchpad<Real> {};
  struct DofsExtScratch : public initializer::VariantScratchpad<Real> {};

  struct FlagScratch : public initializer::Scratchpad<unsigned> {};
  struct QStressNodalScratch : public initializer::VariantScratchpad<Real> {};

  struct ZinvExtra : public initializer::VariantScratchpad<Real> {};

  struct Integrals : public initializer::VariantVariable<DofsArray> {};

  struct LTSVarmap : public initializer::SpecificVarmap<Dofs,
                                                        DofsHalo,
                                                        DofsAne,
                                                        StepIntegrals,
                                                        AccumulatedIntegrals,
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
                                                        Buffers,
                                                        StepIntegralsDevice,
                                                        AccumulatedIntegralsDevice,
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
                                                        QStressNodalScratch,
                                                        Integrals,
                                                        EnergyData,
                                                        ZinvExtra> {};

  using Storage = initializer::Storage<LTSVarmap>;
  using Layer = initializer::Layer<LTSVarmap>;
  template <typename Cfg>
  using Ref = initializer::Layer<LTSVarmap>::CellRef<Cfg>;
  using Backmap = initializer::StorageBackmap<Cell::NumFaces>;

  static void addTo(Storage& storage, const SimulationSettings& settings) {
    using namespace initializer;
    LayerMask plasticityMask;
    if (settings.plasticity) {
      plasticityMask = LayerMask(Ghost);
    } else {
      plasticityMask = LayerMask(Ghost) | LayerMask(Copy) | LayerMask(Interior);
    }
    LayerMask integralMask;
    if (settings.integrate) {
      integralMask = LayerMask(Ghost);
    } else {
      integralMask = LayerMask(Ghost) | LayerMask(Copy) | LayerMask(Interior);
    }

    storage.add<Dofs>(LayerMask(Ghost), PagesizeHeap, allocationModeWP(AllocationPreset::Dofs));
    storage.add<DofsHalo>(LayerMask(Copy) | LayerMask(Interior),
                          PagesizeHeap,
                          allocationModeWP(AllocationPreset::Dofs));

    if (kernels::size<tensor::Qane<Config>>() > 0) {
      storage.add<DofsAne>(
          LayerMask(Ghost), PagesizeHeap, allocationModeWP(AllocationPreset::Dofs));
    } else {
      storage.add<DofsAne>(LayerMask(Ghost) | LayerMask(Copy) | LayerMask(Interior),
                           PagesizeHeap,
                           allocationModeWP(AllocationPreset::Dofs));
    }

    storage.add<StepIntegrals>(
        LayerMask(), Alignment, allocationModeWP(AllocationPreset::TimedofsConstant), true);
    storage.add<AccumulatedIntegrals>(
        LayerMask(), Alignment, allocationModeWP(AllocationPreset::TimedofsConstant), true);
    storage.add<Derivatives>(
        LayerMask(), Alignment, allocationModeWP(AllocationPreset::TimedofsConstant), true);
    storage.add<CellInformation>(
        LayerMask(), Alignment, allocationModeWP(AllocationPreset::Constant), true);
    storage.add<SecondaryInformation>(LayerMask(), Alignment, AllocationMode::HostOnly, true);
    storage.add<FaceNeighbors>(
        LayerMask(Ghost), Alignment, allocationModeWP(AllocationPreset::TimedofsConstant), true);
    storage.add<LocalIntegration>(
        LayerMask(Ghost), Alignment, allocationModeWP(AllocationPreset::ConstantShared), true);
    storage.add<NeighboringIntegration>(
        LayerMask(Ghost), Alignment, allocationModeWP(AllocationPreset::ConstantShared), true);
    storage.add<MaterialData>(LayerMask(), Alignment, AllocationMode::HostOnly, true);
    storage.add<Material>(LayerMask(Ghost), Alignment, AllocationMode::HostOnly, true);
    storage.add<Plasticity>(
        plasticityMask, Alignment, allocationModeWP(AllocationPreset::Plasticity), true);
    storage.add<DRMapping>(
        LayerMask(Ghost), Alignment, allocationModeWP(AllocationPreset::Constant), true);
    storage.add<BoundaryMapping>(
        LayerMask(Ghost), Alignment, allocationModeWP(AllocationPreset::Constant), true);
    storage.add<PStrain>(
        plasticityMask, PagesizeHeap, allocationModeWP(AllocationPreset::PlasticityData));
    storage.add<FaceDisplacements>(LayerMask(Ghost), PagesizeHeap, AllocationMode::HostOnly, true);

    // TODO(David): remove/rename "constant" flag (the data is temporary; and copying it for IO is
    // handled differently)
    storage.add<Buffers>(
        LayerMask(), PagesizeHeap, allocationModeWP(AllocationPreset::Timebucket), true);

    storage.add<StepIntegralsDevice>(LayerMask(), Alignment, AllocationMode::HostOnly, true);
    storage.add<AccumulatedIntegralsDevice>(LayerMask(), Alignment, AllocationMode::HostOnly, true);
    storage.add<DerivativesDevice>(LayerMask(), Alignment, AllocationMode::HostOnly, true);
    storage.add<FaceDisplacementsDevice>(
        LayerMask(Ghost), Alignment, AllocationMode::HostOnly, true);
    storage.add<FaceNeighborsDevice>(LayerMask(Ghost), Alignment, AllocationMode::HostOnly, true);
    storage.add<DRMappingDevice>(LayerMask(Ghost), Alignment, AllocationMode::HostOnly, true);
    storage.add<BoundaryMappingDevice>(LayerMask(Ghost), Alignment, AllocationMode::HostOnly, true);

    storage.add<EnergyData>(LayerMask(Ghost), Alignment, AllocationMode::HostOnly, true);
    storage.add<Integrals>(integralMask, Alignment, allocationModeWP(AllocationPreset::Dofs));

    if constexpr (isDeviceOn()) {
      const auto mode = AllocationMode::DeviceOnly;

      storage.add<DerivativesExtScratch>(LayerMask(), Alignment, mode);
      storage.add<DerivativesAneScratch>(LayerMask(), Alignment, mode);
      storage.add<IDofsAneScratch>(LayerMask(), Alignment, mode);
      storage.add<DofsExtScratch>(LayerMask(), Alignment, mode);
      storage.add<IntegratedDofsScratch>(LayerMask(), Alignment, mode);
      storage.add<DerivativesScratch>(LayerMask(), Alignment, mode);
      storage.add<NodalAvgDisplacements>(LayerMask(), Alignment, mode);
      storage.add<AnalyticScratch>(LayerMask(), Alignment, AllocationMode::HostDevicePinned);

      storage.add<FlagScratch>(LayerMask(), Alignment, mode);
      storage.add<QStressNodalScratch>(LayerMask(), Alignment, mode);

      storage.add<ZinvExtra>(LayerMask(), Alignment, AllocationMode::HostDevicePinned);
    }
  }

  static void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                          Storage& storage) {
    manager.registerData<Dofs>("dofs", storage);
    if constexpr (kernels::size<tensor::Qane<Config>>() > 0) {
      manager.registerData<DofsAne>("dofsAne", storage);
    }
    // check plasticity usage over the layer mask (for now)
    if (storage.info<Plasticity>().mask == initializer::LayerMask(Ghost)) {
      manager.registerData<Plasticity>("pstrain", storage);
    }
  }
};

} // namespace seissol

#endif // SEISSOL_SRC_MEMORY_DESCRIPTOR_LTS_H_
