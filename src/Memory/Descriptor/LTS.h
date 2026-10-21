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

namespace seissol::tensor {
struct Qane;
struct fluxBuffer;
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

  struct Dofs : public initializer::Variable<real[tensor::Q::size()]> {};
  struct DofsHalo : public initializer::Variable<real[tensor::Q::size()]> {};
  // size is zero if Qane is not defined
  struct DofsAne : public initializer::Variable<real[zeroGuard(kernels::size<tensor::Qane>())]> {};
  /// Where the cell writes the trace of its time integral over its step on each of its faces: for
  /// its own corrector, and for the neighbor across the face where that runs at the same time step,
  /// which reads it from there (FaceBuffers points to the same place then). In the Buffers bucket.
  struct FaceTraces : public initializer::Variable<std::array<real*, Cell::NumFaces>> {};
  /// What the cell writes on each face for the neighbor across it (LtsSetup::faceBuffer), null
  /// where it writes nothing there. In the Buffers bucket.
  struct FaceBuffers : public initializer::Variable<std::array<real*, Cell::NumFaces>> {};
  struct CellInformation : public initializer::Variable<CellLocalInformation> {};
  struct SecondaryInformation : public initializer::Variable<SecondaryCellLocalInformation> {};
  /// What the neighbor across each face writes there for this cell: its FaceBuffers of that face.
  struct FaceNeighbors : public initializer::Variable<std::array<real*, Cell::NumFaces>> {};
  struct LocalIntegration : public initializer::Variable<LocalIntegrationData> {};
  struct NeighboringIntegration : public initializer::Variable<NeighboringIntegrationData> {};
  struct MaterialData : public initializer::Variable<model::MaterialT> {};
  /// How many points the material is sampled at inside a cell. Which set that
  /// is, the build decides; it follows the plastic strain unless asked
  /// otherwise.
  static constexpr std::size_t MaterialNodes = tensor::materialNodes::Shape[0];
  struct NodalMaterialData
      : public initializer::Variable<std::array<model::MaterialT, MaterialNodes>> {};
  struct Material : public initializer::Variable<CellMaterialData> {};
  struct Plasticity : public initializer::Variable<seissol::model::PlasticityData> {};
  struct DRMapping : public initializer::Variable<std::array<CellDRMapping, Cell::NumFaces>> {};
  struct BoundaryMapping
      : public initializer::Variable<std::array<CellBoundaryMapping, Cell::NumFaces>> {};
  struct PStrain : public initializer::Variable<
                       real[tensor::QStressNodal::size() + tensor::QEtaNodal::size()]> {};
  struct FaceDisplacements : public initializer::Variable<std::array<real*, Cell::NumFaces>> {};
  struct Buffers : public initializer::Bucket<real> {};
  /// A trace of zeros, which the corrector of a face without a neighbor reads as the neighbor's
  /// (see kernels::applyFaceCorrectors). One per layer.
  struct ZeroFaceTrace : public initializer::Bucket<real> {};

  struct FaceTracesDevice : public initializer::Variable<std::array<real*, Cell::NumFaces>> {};
  struct FaceBuffersDevice : public initializer::Variable<std::array<real*, Cell::NumFaces>> {};
  struct FaceNeighborsDevice : public initializer::Variable<std::array<real*, Cell::NumFaces>> {};
  struct FaceDisplacementsDevice : public initializer::Variable<std::array<real*, Cell::NumFaces>> {
  };
  struct DRMappingDevice : public initializer::Variable<std::array<CellDRMapping, Cell::NumFaces>> {
  };
  struct BoundaryMappingDevice
      : public initializer::Variable<std::array<CellBoundaryMapping, Cell::NumFaces>> {};

  struct EnergyData : public initializer::Variable<model::MaterialT::EnergyData> {};

  struct IntegratedDofsScratch : public initializer::Scratchpad<real> {};
  struct DerivativesScratch : public initializer::Scratchpad<real> {};
  struct NodalAvgDisplacements : public initializer::Scratchpad<real> {};
  struct AnalyticScratch : public initializer::Scratchpad<real> {};
  struct DerivativesExtScratch : public initializer::Scratchpad<real> {};
  struct DerivativesAneScratch : public initializer::Scratchpad<real> {};
  struct IDofsAneScratch : public initializer::Scratchpad<real> {};
  struct DofsExtScratch : public initializer::Scratchpad<real> {};

  struct FlagScratch : public initializer::Scratchpad<unsigned> {};
  struct QStressNodalScratch : public initializer::Scratchpad<real> {};

  struct ZinvExtra : public initializer::Scratchpad<real> {};

  struct Integrals : public initializer::Variable<real[tensor::Q::size()]> {};

  /// What the faces of a weight-adjusted cell lift without material within one step, until the
  /// cell applies the inverse weight of its medium to it (see WeightAdjustable). It has the shape
  /// of what the faces lift into: the field, or the extended field of a solver that keeps the
  /// memory variables apart. Every cell of a rank has one where any of them is weight-adjusted; one
  /// that is not leaves it alone.
  struct FluxBuffer
      : public initializer::Variable<real[zeroGuard(kernels::size<tensor::fluxBuffer>())]> {};

  struct LTSVarmap : public initializer::SpecificVarmap<Dofs,
                                                        DofsHalo,
                                                        DofsAne,
                                                        FaceTraces,
                                                        FaceBuffers,
                                                        CellInformation,
                                                        SecondaryInformation,
                                                        FaceNeighbors,
                                                        LocalIntegration,
                                                        NeighboringIntegration,
                                                        Material,
                                                        MaterialData,
                                                        NodalMaterialData,
                                                        Plasticity,
                                                        DRMapping,
                                                        BoundaryMapping,
                                                        PStrain,
                                                        FaceDisplacements,
                                                        Buffers,
                                                        ZeroFaceTrace,
                                                        FaceTracesDevice,
                                                        FaceBuffersDevice,
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
                                                        FluxBuffer,
                                                        EnergyData,
                                                        ZinvExtra> {};

  using Storage = initializer::Storage<LTSVarmap>;
  using Layer = initializer::Layer<LTSVarmap>;
  using Ref = initializer::Layer<LTSVarmap>::CellRef;
  using Backmap = initializer::StorageBackmap<Cell::NumFaces>;

  static void addTo(Storage& storage, const SimulationSettings& settings) {
    using namespace initializer;
    LayerMask plasticityMask;
    if (settings.plasticity) {
      plasticityMask = LayerMask(Ghost);
    } else {
      plasticityMask = LayerMask(Ghost) | LayerMask(Copy) | LayerMask(Interior);
    }
    LayerMask materialNodalMask;
    if (settings.materialNodal) {
      // the ghost layer carries the samples too: a fault face reads the material
      // of both its cells at its own points, and one of them can sit on another
      // rank
      materialNodalMask = LayerMask();
    } else {
      materialNodalMask = LayerMask(Ghost) | LayerMask(Copy) | LayerMask(Interior);
    }
    LayerMask integralMask;
    if (settings.integrate) {
      integralMask = LayerMask(Ghost);
    } else {
      integralMask = LayerMask(Ghost) | LayerMask(Copy) | LayerMask(Interior);
    }
    LayerMask fluxBufferMask;
    if (settings.weightAdjusted) {
      fluxBufferMask = LayerMask(Ghost);
    } else {
      fluxBufferMask = LayerMask(Ghost) | LayerMask(Copy) | LayerMask(Interior);
    }

    storage.add<Dofs>(LayerMask(Ghost), PagesizeHeap, allocationModeWP(AllocationPreset::Dofs));
    storage.add<DofsHalo>(LayerMask(Copy) | LayerMask(Interior),
                          PagesizeHeap,
                          allocationModeWP(AllocationPreset::Dofs));

    if (kernels::size<tensor::Qane>() > 0) {
      storage.add<DofsAne>(
          LayerMask(Ghost), PagesizeHeap, allocationModeWP(AllocationPreset::Dofs));
    } else {
      storage.add<DofsAne>(LayerMask(Ghost) | LayerMask(Copy) | LayerMask(Interior),
                           PagesizeHeap,
                           allocationModeWP(AllocationPreset::Dofs));
    }

    storage.add<FaceTraces>(
        LayerMask(), Alignment, allocationModeWP(AllocationPreset::TimedofsConstant), true);
    storage.add<FaceBuffers>(
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
    storage.add<NodalMaterialData>(materialNodalMask, Alignment, AllocationMode::HostOnly, true);
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
    storage.add<ZeroFaceTrace>(
        LayerMask(Ghost), Alignment, allocationModeWP(AllocationPreset::Timebucket), true);

    storage.add<FaceTracesDevice>(LayerMask(), Alignment, AllocationMode::HostOnly, true);
    storage.add<FaceBuffersDevice>(LayerMask(), Alignment, AllocationMode::HostOnly, true);
    storage.add<FaceDisplacementsDevice>(
        LayerMask(Ghost), Alignment, AllocationMode::HostOnly, true);
    storage.add<FaceNeighborsDevice>(LayerMask(Ghost), Alignment, AllocationMode::HostOnly, true);
    storage.add<DRMappingDevice>(LayerMask(Ghost), Alignment, AllocationMode::HostOnly, true);
    storage.add<BoundaryMappingDevice>(LayerMask(Ghost), Alignment, AllocationMode::HostOnly, true);

    storage.add<EnergyData>(LayerMask(Ghost), Alignment, AllocationMode::HostOnly, true);
    storage.add<Integrals>(integralMask, Alignment, allocationModeWP(AllocationPreset::Dofs));
    storage.add<FluxBuffer>(fluxBufferMask, PagesizeHeap, allocationModeWP(AllocationPreset::Dofs));

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
    if constexpr (kernels::size<tensor::Qane>() > 0) {
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
