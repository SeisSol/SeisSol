// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_MEMORY_DESCRIPTOR_LTS_H_
#define SEISSOL_SRC_MEMORY_DESCRIPTOR_LTS_H_

#include "IO/Instance/Checkpoint/CheckpointManager.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"
#include "Model/Plasticity.h"
#include "generated_code/tensor.h"
#include <Initializer/CellLocalInformation.h>
#include <Parallel/Helper.h>

#ifdef ACL_DEVICE
#include "Parallel/Helper.h"
#endif

namespace seissol::tensor {
class Qane;
} // namespace seissol::tensor

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
  struct Dofs : public initializer::Variable<real[tensor::Q::size()]> {};
  // size is zero if Qane is not defined
  struct DofsAne
      : public initializer::Variable<real[zeroLengthArrayHandler(kernels::size<tensor::Qane>())]> {
  };
  struct Buffers : public initializer::Variable<real*> {};
  struct Derivatives : public initializer::Variable<real*> {};
  struct CellInformation : public initializer::Variable<CellLocalInformation> {};
  struct SecondaryInformation : public initializer::Variable<SecondaryCellLocalInformation> {};
  struct FaceNeighbors : public initializer::Variable<real* [Cell::NumFaces]> {};
  struct LocalIntegration : public initializer::Variable<LocalIntegrationData> {};
  struct NeighboringIntegration : public initializer::Variable<NeighboringIntegrationData> {};
  struct Material : public initializer::Variable<CellMaterialData> {};
  struct Plasticity : public initializer::Variable<seissol::model::PlasticityData> {};
  struct DRMapping : public initializer::Variable<CellDRMapping[Cell::NumFaces]> {};
  struct BoundaryMapping : public initializer::Variable<CellBoundaryMapping[Cell::NumFaces]> {};
  struct PStrain
      : public initializer::Variable<real[tensor::QStress::size() + tensor::QEtaModal::size()]> {};
  struct FaceDisplacements : public initializer::Variable<real* [Cell::NumFaces]> {};
  struct BuffersDerivatives : public initializer::Bucket<real> {};

  struct BuffersDevice : public initializer::Variable<real*> {};
  struct DerivativesDevice : public initializer::Variable<real*> {};
  struct FaceNeighborsDevice : public initializer::Variable<real* [Cell::NumFaces]> {};
  struct FaceDisplacementsDevice : public initializer::Variable<real* [Cell::NumFaces]> {};
  struct DRMappingDevice : public initializer::Variable<CellDRMapping[Cell::NumFaces]> {};
  struct BoundaryMappingDevice : public initializer::Variable<CellBoundaryMapping[Cell::NumFaces]> {
  };

  struct IntegratedDofsScratch : public initializer::Scratchpad<real> {};
  struct DerivativesScratch : public initializer::Scratchpad<real> {};
  struct NodalAvgDisplacements : public initializer::Scratchpad<real> {};
  struct AnalyticScratch : public initializer::Scratchpad<real> {};
  struct DerivativesExtScratch : public initializer::Scratchpad<real> {};
  struct DerivativesAneScratch : public initializer::Scratchpad<real> {};
  struct IDofsAneScratch : public initializer::Scratchpad<real> {};
  struct DofsExtScratch : public initializer::Scratchpad<real> {};

  struct FlagScratch : public initializer::Scratchpad<unsigned> {};
  struct PrevDofsScratch : public initializer::Scratchpad<real> {};
  struct QEtaNodalScratch : public initializer::Scratchpad<real> {};
  struct QStressNodalScratch : public initializer::Scratchpad<real> {};

  struct RotateDisplacementToFaceNormalScratch : public initializer::Scratchpad<real> {};
  struct RotateDisplacementToGlobalScratch : public initializer::Scratchpad<real> {};
  struct RotatedFaceDisplacementScratch : public initializer::Scratchpad<real> {};
  struct DofsFaceNodalScratch : public initializer::Scratchpad<real> {};
  struct PrevCoefficientsScratch : public initializer::Scratchpad<real> {};
  struct DofsFaceBoundaryNodalScratch : public initializer::Scratchpad<real> {};

  struct Integrals : public initializer::Variable<real> {};

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

  using Tree = initializer::LTSTree<LTSVarmap>;
  using Layer = initializer::Layer<LTSVarmap>;
  using Ref = initializer::Layer<LTSVarmap>::CellRef;

  static void addTo(Tree& tree, bool usePlasticity) {
    using namespace initializer;
    LayerMask plasticityMask;
    if (usePlasticity) {
      plasticityMask = LayerMask(Ghost);
    } else {
      plasticityMask = LayerMask(Ghost) | LayerMask(Copy) | LayerMask(Interior);
    }

    tree.add<Dofs>(LayerMask(Ghost), PagesizeHeap, allocationModeWP(AllocationPreset::Dofs));
    if (kernels::size<tensor::Qane>() > 0) {
      tree.add<DofsAne>(LayerMask(Ghost), PagesizeHeap, allocationModeWP(AllocationPreset::Dofs));
    } else {
      tree.add<DofsAne>(LayerMask(Ghost) | LayerMask(Copy) | LayerMask(Interior),
                        PagesizeHeap,
                        allocationModeWP(AllocationPreset::Dofs));
    }
    tree.add<Buffers>(LayerMask(), 1, allocationModeWP(AllocationPreset::TimedofsConstant), true);
    tree.add<Derivatives>(
        LayerMask(), 1, allocationModeWP(AllocationPreset::TimedofsConstant), true);
    tree.add<CellInformation>(LayerMask(), 1, allocationModeWP(AllocationPreset::Constant), true);
    tree.add<SecondaryInformation>(LayerMask(), 1, AllocationMode::HostOnly, true);
    tree.add<FaceNeighbors>(
        LayerMask(Ghost), 1, allocationModeWP(AllocationPreset::TimedofsConstant), true);
    tree.add<LocalIntegration>(
        LayerMask(Ghost), 1, allocationModeWP(AllocationPreset::ConstantShared), true);
    tree.add<NeighboringIntegration>(
        LayerMask(Ghost), 1, allocationModeWP(AllocationPreset::ConstantShared), true);
    tree.add<Material>(LayerMask(Ghost), 1, AllocationMode::HostOnly, true);
    tree.add<Plasticity>(plasticityMask, 1, allocationModeWP(AllocationPreset::Plasticity), true);
    tree.add<DRMapping>(LayerMask(Ghost), 1, allocationModeWP(AllocationPreset::Constant), true);
    tree.add<BoundaryMapping>(
        LayerMask(Ghost), 1, allocationModeWP(AllocationPreset::Constant), true);
    tree.add<PStrain>(
        plasticityMask, PagesizeHeap, allocationModeWP(AllocationPreset::PlasticityData));
    tree.add<FaceDisplacements>(LayerMask(Ghost), PagesizeHeap, AllocationMode::HostOnly, true);

    // TODO(David): remove/rename "constant" flag (the data is temporary; and copying it for IO is
    // handled differently)
    tree.add<BuffersDerivatives>(
        LayerMask(), PagesizeHeap, allocationModeWP(AllocationPreset::Timebucket), true);

    tree.add<BuffersDevice>(LayerMask(), 1, AllocationMode::HostOnly, true);
    tree.add<DerivativesDevice>(LayerMask(), 1, AllocationMode::HostOnly, true);
    tree.add<FaceDisplacementsDevice>(LayerMask(Ghost), 1, AllocationMode::HostOnly, true);
    tree.add<FaceNeighborsDevice>(LayerMask(Ghost), 1, AllocationMode::HostOnly, true);
    tree.add<DRMappingDevice>(LayerMask(Ghost), 1, AllocationMode::HostOnly, true);
    tree.add<BoundaryMappingDevice>(LayerMask(Ghost), 1, AllocationMode::HostOnly, true);

    if constexpr (isDeviceOn()) {
      tree.add<DerivativesExtScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add<DerivativesAneScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add<IDofsAneScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add<DofsExtScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add<IntegratedDofsScratch>(LayerMask(), 1, AllocationMode::HostDeviceSplit);
      tree.add<DerivativesScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add<NodalAvgDisplacements>(LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add<AnalyticScratch>(LayerMask(), 1, AllocationMode::HostDevicePinned);

      tree.add<FlagScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add<PrevDofsScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add<QEtaNodalScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add<QStressNodalScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);

      tree.add<RotateDisplacementToFaceNormalScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add<RotateDisplacementToGlobalScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add<RotatedFaceDisplacementScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add<DofsFaceNodalScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add<PrevCoefficientsScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add<DofsFaceBoundaryNodalScratch>(LayerMask(), 1, AllocationMode::DeviceOnly);
    }
  }

  static void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                          Tree* tree) {
    manager.registerData<Dofs>("dofs", tree);
    if constexpr (kernels::size<tensor::Qane>() > 0) {
      manager.registerData<DofsAne>("dofsAne", tree);
    }
    // check plasticity usage over the layer mask (for now)
    if (tree->info<Plasticity>().mask == initializer::LayerMask(Ghost)) {
      manager.registerData<Plasticity>("pstrain", tree);
    }
  }
};

} // namespace seissol

#endif // SEISSOL_SRC_MEMORY_DESCRIPTOR_LTS_H_
