// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#include "Neighbor.h"

#include "Common/Marker.h"
#include "GeneratedCode/init.h"
#include "Kernels/StarOperands.h"
#include "Model/OperatorLayout.h"
#include "Monitoring/Metric.h"

#include <cassert>
#include <cstddef>
#include <cstring>
#include <iterator>
#include <stdint.h>

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#endif

namespace seissol::kernels::solver::linearckanelastic {

// The neighbouring flux family is indexed by the neighbouring side and the own face. The face
// orientation index is not part of it, since the canonical vertex numbering pins it to zero on
// every interior face.
static_assert(std::size(seissol::kernel::neighborFluxExt::ExecutePtrs) ==
              Cell::NumFaces * Cell::NumFaces);

void Neighbor::setGlobalData(const CompoundGlobalData& global) {
  nfKrnlPrototype_.bindGlobals(*global.onHost);
  drKrnlPrototype_.bindGlobals(*global.onHost);

#ifdef ACL_DEVICE
  deviceNfKrnlPrototype_.bindGlobals(*global.onDevice);
  deviceDrKrnlPrototype_.bindGlobals(*global.onDevice);
#endif
}

void Neighbor::computeNeighborsIntegral(
    LTS::Ref& data,
    const std::array<real*, Cell::NumFaces>& timeIntegrated,
    const std::array<real*, Cell::NumFaces>& faceNeighborsPrefetch) {
#ifndef NDEBUG
  for (std::size_t neighbor = 0; neighbor < Cell::NumFaces; ++neighbor) {
    // alignment of the time integrated dofs (only for linear interior)
    if (data.get<LTS::CellInformation>().faceTypes[neighbor] == FaceType::Regular) {
      assert((reinterpret_cast<uintptr_t>(timeIntegrated[neighbor])) % Vectorsize == 0);
    }
  }
#endif

  const auto& cellDrMapping = data.get<LTS::DRMapping>();

  // alignment of the degrees of freedom
  assert((reinterpret_cast<uintptr_t>(data.get<LTS::Dofs>())) % Vectorsize == 0);

  alignas(PagesizeStack) real Qext[tensor::Qext::size()] = {};

  kernel::neighborFluxExt nfKrnl = nfKrnlPrototype_;
  nfKrnl.Qext = Qext;

  // iterate over faces
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // neighboring cell contribution only for interior faces
    if (data.get<LTS::CellInformation>().faceTypes[face] == FaceType::Regular) {
      assert(data.get<LTS::CellInformation>().faceRelations[face][0] < Cell::NumFaces &&
             data.get<LTS::CellInformation>().faceRelations[face][1] == 0);

      nfKrnl.I = timeIntegrated[face];
      kernels::bindNeighborFluxOperands(
          nfKrnl, data.get<LTS::LocalIntegration>(), data.get<LTS::NeighboringIntegration>(), face);
      nfKrnl._prefetch.I = faceNeighborsPrefetch[face];
      nfKrnl.execute(data.get<LTS::CellInformation>().faceRelations[face][0], face);
    } else if (data.get<LTS::CellInformation>().faceTypes[face] == FaceType::DynamicRupture) {
      assert((reinterpret_cast<uintptr_t>(cellDrMapping[face].godunov)) % Vectorsize == 0);

      dynamicRupture::kernel::nodalFlux drKrnl = drKrnlPrototype_;
      kernels::bindFaultFluxOperands(drKrnl, cellDrMapping[face].fluxSolver);
      drKrnl.QInterpolated = cellDrMapping[face].godunov;
      drKrnl.Qext = Qext;
      drKrnl._prefetch.I = faceNeighborsPrefetch[face];
      drKrnl.execute(cellDrMapping[face].side, cellDrMapping[face].faceRelation);
    }
  }

  kernel::neighbor nKrnl = nKrnlPrototype_;
  nKrnl.Qext = Qext;
  nKrnl.Q = data.get<LTS::Dofs>();
  nKrnl.Qane = data.get<LTS::DofsAne>();
  nKrnl.w = data.get<LTS::NeighboringIntegration>().specific.w;

  nKrnl.execute();
}

std::pair<PerformanceEstimate, PerformanceEstimate>
    Neighbor::metrics(const std::array<FaceType, Cell::NumFaces>& faceTypes,
                      const std::array<std::array<uint8_t, 2>, Cell::NumFaces>& neighboringIndices,
                      const std::array<CellDRMapping, Cell::NumFaces>& cellDrMapping) const {

  PerformanceEstimate regular;
  PerformanceEstimate dr;

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // neighboring cell contribution only for interior faces
    if (faceTypes[face] == FaceType::Regular) {
      assert(neighboringIndices[face][0] < Cell::NumFaces && neighboringIndices[face][1] == 0);

      regular += PerformanceEstimate::fromKernel<seissol::kernel::neighborFluxExt>(
          neighboringIndices[face][0], face);
    } else if (faceTypes[face] == FaceType::DynamicRupture) {
      dr += PerformanceEstimate::fromKernel<dynamicRupture::kernel::nodalFlux>(
          cellDrMapping[face].side, cellDrMapping[face].faceRelation);
    }
  }

  regular += PerformanceEstimate::fromKernel<kernel::neighbor>();

  // legacy memory estimate
  std::uint64_t reals = 0;

  // 4 * tElasticDOFS load, DOFs load, DOFs write
  reals += 4 * tensor::I::size() + 2 * tensor::Q::size() + 2 * tensor::Qane::size();
  // flux solvers load
  reals += 4 * tensor::AminusT::size() + tensor::w::size();

  regular.bytes = reals * sizeof(real);

  return {regular, dr};
}

void Neighbor::computeBatchedNeighborsIntegral(
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& table,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE

  using namespace seissol::recording;
  kernel::gpu_neighborFluxExt neighFluxKrnl = deviceNfKrnlPrototype_;
  dynamicRupture::kernel::gpu_nodalFlux drKrnl = deviceDrKrnlPrototype_;

  {
    ConditionalKey key(KernelNames::Time || KernelNames::Volume);
    if (table.find(key) != table.end()) {
      auto& entry = table[key];
      device.algorithms.setToValue((entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr(),
                                   static_cast<real>(0.0),
                                   tensor::Qext::Size,
                                   (entry.get(inner_keys::Wp::Id::DofsExt))->getSize(),
                                   runtime.stream());
    }
  }

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    runtime.envMany(
        (*FaceRelations::Count) + (*DrFaceRelations::Count), [&](void* stream, size_t i) {
          // regular and periodic
          if (i < (*FaceRelations::Count)) {
            // regular and periodic
            const auto faceRelation = i;

            ConditionalKey key(*KernelNames::NeighborFlux, *FaceKinds::Regular, face, faceRelation);

            if (table.find(key) != table.end()) {
              auto& entry = table[key];

              const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
              neighFluxKrnl.numElements = numElements;

              neighFluxKrnl.Qext = (entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr();
              neighFluxKrnl.I = const_cast<const real**>(
                  (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());
              // the cell's own data is recorded only where the flux reads the
              // rotation of the face from it
              const real** localIntegrationPtrs = nullptr;
              if constexpr (NodalFlux) {
                localIntegrationPtrs = const_cast<const real**>(
                    entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
              }
              kernels::bindNeighborFluxOperandsBatched(
                  neighFluxKrnl,
                  localIntegrationPtrs,
                  const_cast<const real**>(
                      entry.get(inner_keys::Wp::Id::NeighborIntegrationData)->getDeviceDataPtr()),
                  face);

              neighFluxKrnl.streamPtr = stream;
              (neighFluxKrnl.*seissol::kernel::gpu_neighborFluxExt::ExecutePtrs[faceRelation])();
            }
          } else {
            // Dynamic Rupture
            const auto faceRelation = i - (*FaceRelations::Count);

            ConditionalKey key(
                *KernelNames::NeighborFlux, *FaceKinds::DynamicRupture, face, faceRelation);

            if (table.find(key) != table.end()) {
              auto& entry = table[key];

              const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
              drKrnl.numElements = numElements;

              kernels::bindFaultFluxOperandsBatched(
                  drKrnl,
                  const_cast<const real**>(
                      (entry.get(inner_keys::Wp::Id::FluxSolver))->getDeviceDataPtr()));
              drKrnl.QInterpolated = const_cast<const real**>(
                  (entry.get(inner_keys::Wp::Id::Godunov))->getDeviceDataPtr());
              drKrnl.Qext = (entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr();

              // the lift keeps a temporary where the face carries it per point
              real* tmpMem = nullptr;
              if constexpr (seissol::dynamicRupture::kernel::gpu_nodalFlux::
                                TmpMaxMemRequiredInBytes > 0) {
                tmpMem = reinterpret_cast<real*>(device.api->allocMemAsync(
                    seissol::dynamicRupture::kernel::gpu_nodalFlux::TmpMaxMemRequiredInBytes *
                        numElements,
                    stream));
                drKrnl.linearAllocator.initialize(tmpMem);
              }

              drKrnl.streamPtr = stream;
              (drKrnl.*seissol::dynamicRupture::kernel::gpu_nodalFlux::ExecutePtrs[faceRelation])();
              if (tmpMem != nullptr) {
                device.api->freeMemAsync(reinterpret_cast<void*>(tmpMem), stream);
              }
            }
          }
        });
  }

  ConditionalKey key(KernelNames::Time || KernelNames::Volume);
  if (table.find(key) != table.end()) {
    auto& entry = table[key];
    kernel::gpu_neighbor nKrnl = deviceNKrnlPrototype_;
    nKrnl.numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
    nKrnl.Qext =
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr());
    nKrnl.Q = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
    nKrnl.Qane = (entry.get(inner_keys::Wp::Id::DofsAne))->getDeviceDataPtr();
    nKrnl.w = const_cast<const real**>(
        entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
    nKrnl.extraOffset_w = SEISSOL_OFFSET(LocalIntegrationData, specific.w);

    SEISSOL_OFFSET_ASSERT(LocalIntegrationData, specific.w);

    nKrnl.streamPtr = runtime.stream();

    nKrnl.execute();
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

} // namespace seissol::kernels::solver::linearckanelastic
