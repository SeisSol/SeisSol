// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_WAVEPROPAGATION_SOLVER_STORAGE_H_
#define SEISSOL_SRC_WAVEPROPAGATION_SOLVER_STORAGE_H_

#include "Alignment.h"
#include "Core/Common/BasicTypedefs.h"
#include "Core/Common/Constants.h"
#include "Core/Common/Real.h"
#include "Core/Common/Typedefs.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/SimulationSettings.h"
#include "Core/Memory/Tree/LTSTree.h"
#include "Core/Memory/Tree/Layer.h"
#include "GeneratedCode/tensor.h"
#include "WavePropagation/Physics/Datastructures.h"
#include "WavePropagation/Solver/CellData.h"
#include "WavePropagation/Solver/LinearCKAnelastic/Storage.h" // IWYU pragma: export
#include "WavePropagation/Solver/STP/Storage.h"               // IWYU pragma: export

#include <array>
#include <cstdint>
#include <functional>
#include <utility>

/**
 * The variables of the cells of the wave propagation in the storage of the cells (see LTS): their
 * unknowns, their buffers and derivatives, the data of their integration, their materials, and the
 * mappings of their faces. The solvers add what they keep beyond that (see
 * Solver/<solver>/Storage.h).
 */
namespace seissol::LTS {

// The unknowns of a cell and what is derived from them are held in the reals and the layout of
// the configuration of their layer.
template <typename Cfg>
using DofsArray = Real<Cfg>[tensor::Q<Cfg>::size()];
template <typename Cfg>
using RealPtr = Real<Cfg>*;
template <typename Cfg>
using FaceRealPtrs = std::array<Real<Cfg>*, Cell::NumFaces>;
// The data of the integration, the material and the mappings of the faces of a cell are those
// of the configuration of its layer as well.
template <typename Cfg>
using FaceDRMappings = std::array<CellDRMapping<Cfg>, Cell::NumFaces>;
template <typename Cfg>
using FaceBoundaryMappings = std::array<CellBoundaryMapping<Cfg>, Cell::NumFaces>;
template <typename Cfg>
using EnergyDataOf = model::MaterialOf<Cfg>::template EnergyData<Cfg>;

struct Dofs : public initializer::VariantVariable<DofsArray> {};
struct DofsHalo : public initializer::VariantVariable<DofsArray> {};
struct StepIntegrals : public initializer::VariantVariable<RealPtr> {};
struct AccumulatedIntegrals : public initializer::VariantVariable<RealPtr> {};
struct Derivatives : public initializer::VariantVariable<RealPtr> {};
// The buffers or derivatives of the neighbors, which hold them in the reals of their own
// configuration.
struct FaceNeighbors : public initializer::Variable<std::array<void*, Cell::NumFaces>> {};
struct LocalIntegration : public initializer::VariantVariable<LocalIntegrationData> {};
struct NeighboringIntegration : public initializer::VariantVariable<NeighboringIntegrationData> {};
// per face towards a neighbor of another family in another configuration: the coupling matrix
// of the face, the contribution of the neighbor to the Godunov state in the quantities of the
// cell (see kernels::ConfigBoundary), in the bucket CouplingMatrices; null for the other faces
struct Coupling : public initializer::VariantVariable<FaceRealPtrs> {};
struct CouplingMatrices : public initializer::VariantBucket<Real> {};
struct MaterialData : public initializer::VariantVariable<model::MaterialOf> {};
struct Material : public initializer::Variable<CellMaterialData> {};
struct DRMapping : public initializer::VariantVariable<FaceDRMappings> {};
struct BoundaryMapping : public initializer::VariantVariable<FaceBoundaryMappings> {};
struct FaceDisplacements : public initializer::VariantVariable<FaceRealPtrs> {};
// The buffers and derivatives of the cells (see StepIntegrals, AccumulatedIntegrals and
// Derivatives): of the interior in Buffers; of the copy and ghost layers, which MPI sends and
// receives, in BuffersHalo, which SEISSOL_USM_MPI makes unified memory (see
// AllocationPreset::TimebucketHalo). withBuffersOf() picks the one of a layer.
struct Buffers : public initializer::VariantBucket<Real> {};
struct BuffersHalo : public initializer::VariantBucket<Real> {};

struct StepIntegralsDevice : public initializer::VariantVariable<RealPtr> {};
struct AccumulatedIntegralsDevice : public initializer::VariantVariable<RealPtr> {};
struct DerivativesDevice : public initializer::VariantVariable<RealPtr> {};
struct FaceNeighborsDevice : public initializer::Variable<std::array<void*, Cell::NumFaces>> {};
struct FaceDisplacementsDevice : public initializer::VariantVariable<FaceRealPtrs> {};
struct CouplingDevice : public initializer::VariantVariable<FaceRealPtrs> {};
struct DRMappingDevice : public initializer::VariantVariable<FaceDRMappings> {};
struct BoundaryMappingDevice : public initializer::VariantVariable<FaceBoundaryMappings> {};

struct EnergyData : public initializer::VariantVariable<EnergyDataOf> {};

struct IntegratedDofsScratch : public initializer::VariantScratchpad<Real> {};
struct DerivativesScratch : public initializer::VariantScratchpad<Real> {};
struct NodalAvgDisplacements : public initializer::VariantScratchpad<Real> {};
struct AnalyticScratch : public initializer::VariantScratchpad<Real> {};

// the time integrals of the neighbors of other configurations, in their configuration, in the
// canonical form, and converted (see kernels::ConfigBoundary); in bytes, as they mix reals
struct ConfigBoundaryScratch : public initializer::Scratchpad<std::uint8_t> {};

struct Integrals : public initializer::VariantVariable<DofsArray> {};

/// Calls `function` with the bucket that holds the buffers and derivatives of the cells of
/// `layer`: Buffers in the interior, BuffersHalo in the copy and ghost layers.
template <typename F>
void withBuffersOf(const Layer& layer, F&& function) {
  if (layer.getIdentifier().halo == HaloType::Interior) {
    std::invoke(std::forward<F>(function), Buffers{});
  } else {
    std::invoke(std::forward<F>(function), BuffersHalo{});
  }
}

/// Adds the variables of the cells of the wave propagation to `storage`, before it gets
/// allocated; with those of the core (see addTo()).
inline void addWavePropagation(Storage& storage, const SimulationSettings& settings) {
  using namespace initializer;
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

  // what the solvers keep for themselves
  kernels::solver::linearckanelastic::addStorage(storage);
  kernels::solver::stp::addStorage(storage);

  storage.add<StepIntegrals>(
      LayerMask(), Alignment, allocationModeWP(AllocationPreset::TimedofsConstant), true);
  storage.add<AccumulatedIntegrals>(
      LayerMask(), Alignment, allocationModeWP(AllocationPreset::TimedofsConstant), true);
  storage.add<Derivatives>(
      LayerMask(), Alignment, allocationModeWP(AllocationPreset::TimedofsConstant), true);
  storage.add<FaceNeighbors>(
      LayerMask(Ghost), Alignment, allocationModeWP(AllocationPreset::TimedofsConstant), true);
  storage.add<LocalIntegration>(
      LayerMask(Ghost), Alignment, allocationModeWP(AllocationPreset::ConstantShared), true);
  storage.add<NeighboringIntegration>(
      LayerMask(Ghost), Alignment, allocationModeWP(AllocationPreset::ConstantShared), true);
  storage.add<Coupling>(LayerMask(Ghost), Alignment, AllocationMode::HostOnly, true);
  storage.add<CouplingMatrices>(
      LayerMask(Ghost), Alignment, allocationModeWP(AllocationPreset::ConstantShared), true);
  storage.add<MaterialData>(LayerMask(), Alignment, AllocationMode::HostOnly, true);
  storage.add<Material>(LayerMask(Ghost), Alignment, AllocationMode::HostOnly, true);
  storage.add<DRMapping>(
      LayerMask(Ghost), Alignment, allocationModeWP(AllocationPreset::Constant), true);
  storage.add<BoundaryMapping>(
      LayerMask(Ghost), Alignment, allocationModeWP(AllocationPreset::Constant), true);
  storage.add<FaceDisplacements>(LayerMask(Ghost), PagesizeHeap, AllocationMode::HostOnly, true);

  // TODO(David): remove/rename "constant" flag (the data is temporary; and copying it for IO is
  // handled differently)
  storage.add<Buffers>(LayerMask(Ghost) | LayerMask(Copy),
                       PagesizeHeap,
                       allocationModeWP(AllocationPreset::Timebucket),
                       true);
  storage.add<BuffersHalo>(
      LayerMask(Interior), PagesizeHeap, allocationModeWP(AllocationPreset::TimebucketHalo), true);

  storage.add<StepIntegralsDevice>(LayerMask(), Alignment, AllocationMode::HostOnly, true);
  storage.add<AccumulatedIntegralsDevice>(LayerMask(), Alignment, AllocationMode::HostOnly, true);
  storage.add<DerivativesDevice>(LayerMask(), Alignment, AllocationMode::HostOnly, true);
  storage.add<FaceDisplacementsDevice>(LayerMask(Ghost), Alignment, AllocationMode::HostOnly, true);
  storage.add<CouplingDevice>(LayerMask(Ghost), Alignment, AllocationMode::HostOnly, true);
  storage.add<FaceNeighborsDevice>(LayerMask(Ghost), Alignment, AllocationMode::HostOnly, true);
  storage.add<DRMappingDevice>(LayerMask(Ghost), Alignment, AllocationMode::HostOnly, true);
  storage.add<BoundaryMappingDevice>(LayerMask(Ghost), Alignment, AllocationMode::HostOnly, true);

  storage.add<EnergyData>(LayerMask(Ghost), Alignment, AllocationMode::HostOnly, true);
  storage.add<Integrals>(integralMask, Alignment, allocationModeWP(AllocationPreset::Dofs));

  if constexpr (isDeviceOn()) {
    const auto mode = AllocationMode::DeviceOnly;

    storage.add<IntegratedDofsScratch>(LayerMask(), Alignment, mode);
    storage.add<DerivativesScratch>(LayerMask(), Alignment, mode);
    storage.add<NodalAvgDisplacements>(LayerMask(), Alignment, mode);
    storage.add<AnalyticScratch>(LayerMask(), Alignment, AllocationMode::HostDevicePinned);
    storage.add<ConfigBoundaryScratch>(LayerMask(), Alignment, mode);
  }
}

/// Sets the sizes of the scratchpads that the solver of the configuration `Cfg` keeps for itself,
/// in the layer `layer` of its cells.
template <typename Cfg>
void setSolverScratchpadSizes(Layer& layer) {
  if constexpr (Cfg::Solver == SolverType::LinearCKAnelastic) {
    kernels::solver::linearckanelastic::setScratchpadSizes<Cfg>(layer);
  } else if constexpr (Cfg::Solver == SolverType::STP) {
    kernels::solver::stp::setScratchpadSizes<Cfg>(layer);
  }
}

} // namespace seissol::LTS

#endif // SEISSOL_SRC_WAVEPROPAGATION_SOLVER_STORAGE_H_
