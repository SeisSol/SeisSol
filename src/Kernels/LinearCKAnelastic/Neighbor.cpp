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
#include "Kernels/FaceTraces.h"
#include "Kernels/StarOperands.h"
#include "Memory/Descriptor/LTS.h"
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

// The corrector family is indexed by the face and by the side of it the cell is on, which says in
// whose parametrisation the traces of the face are written. The face orientation index is not part
// of it, since the canonical vertex numbering pins it to zero on every interior face.
static_assert(std::size(seissol::kernel::faceCorrector::ExecutePtrs) == Cell::NumFaces * 2);

void Neighbor::setGlobalData(const CompoundGlobalData& global) {
  correctorKrnlPrototype_.bindGlobals(*global.onHost);
  drKrnlPrototype_.bindGlobals(*global.onHost);
  weightAdjustment_.setGlobalData(global);

#ifdef ACL_DEVICE
  deviceCorrectorsKrnlPrototype_.bindGlobals(*global.onDevice);
  deviceDrKrnlPrototype_.bindGlobals(*global.onDevice);
#endif
}

void Neighbor::computeNeighborsIntegral(
    LTS::Ref& data, const std::array<const real*, Cell::NumFaces>& neighborTraces) {
  const auto& cellDrMapping = data.get<LTS::DRMapping>();

  // alignment of the degrees of freedom
  assert((reinterpret_cast<uintptr_t>(data.get<LTS::Dofs>())) % Vectorsize == 0);

  alignas(PagesizeStack) real Qext[tensor::Qext::size()] = {};

  // a weight-adjusted cell collects what its faces lift in its flux buffer, which the local
  // integral started, and adds the inverse weight applied to it to the extended field once they
  // are all in (WeightAdjustable)
  const bool weightAdjusted = WeightAdjustable && data.get<LTS::CellInformation>().weightAdjusted;
  real* fluxTarget = weightAdjusted ? data.get<LTS::FluxBuffer>() : Qext;

  kernel::faceCorrector correctorKrnl = correctorKrnlPrototype_;
  correctorKrnl.Qext = fluxTarget;
  kernels::applyFaceCorrectors(correctorKrnl, neighborTraces, data);

  // the flux the fault supplies on its faces
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    if (data.get<LTS::CellInformation>().faceTypes[face] == FaceType::DynamicRupture) {
      assert((reinterpret_cast<uintptr_t>(cellDrMapping[face].godunov)) % Vectorsize == 0);

      dynamicRupture::kernel::nodalFlux drKrnl = drKrnlPrototype_;
      kernels::bindFaultFluxOperands(drKrnl, cellDrMapping[face].fluxSolver);
      drKrnl.QInterpolated = cellDrMapping[face].godunov;
      drKrnl.Qext = fluxTarget;
      drKrnl.execute(cellDrMapping[face].side, cellDrMapping[face].faceRelation);
    }
  }

  if (weightAdjusted) {
    weightAdjustment_.apply(data, [&](auto& krnl) { krnl.Qext = Qext; });
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
                      const std::array<CellDRMapping, Cell::NumFaces>& cellDrMapping) const {

  PerformanceEstimate regular;
  PerformanceEstimate dr;

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // the corrector of the face; both sides of a face cost the same
    if (appliesFaceCorrector(faceTypes[face])) {
      regular += PerformanceEstimate::fromKernel<seissol::kernel::faceCorrector>(face, 0);
    }
    if (faceTypes[face] == FaceType::DynamicRupture) {
      dr += PerformanceEstimate::fromKernel<dynamicRupture::kernel::nodalFlux>(
          cellDrMapping[face].side, cellDrMapping[face].faceRelation);
    }
  }

  regular += PerformanceEstimate::fromKernel<kernel::neighbor>();

  // legacy memory estimate
  std::uint64_t reals = 0;

  // 8 traces load, DOFs load, DOFs write
  reals += 2 * Cell::NumFaces * tensor::faceTrace::size() + 2 * tensor::Q::size() +
           2 * tensor::Qane::size();
  // relaxation load
  reals += tensor::w::size();

  regular.bytes = reals * sizeof(real);

  return {regular, dr};
}

void Neighbor::computeBatchedNeighborsIntegral(
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& table,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE

  using namespace seissol::recording;

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

  // the corrector of every face, into what the faces of the cells lift into (FluxTarget)
  kernels::applyFaceCorrectorsBatched(
      deviceCorrectorsKrnlPrototype_, table, runtime, [](auto& krnl, auto& entry) {
        krnl.Qext = (entry.get(inner_keys::Wp::Id::FluxTarget))->getDeviceDataPtr();
      });

  // the flux the fault supplies on its faces; the face relations of one face are different cells
  dynamicRupture::kernel::gpu_nodalFlux drKrnl = deviceDrKrnlPrototype_;
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    runtime.envMany((*DrFaceRelations::Count), [&](void* stream, size_t faceRelation) {
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
        drKrnl.QInterpolated =
            const_cast<const real**>((entry.get(inner_keys::Wp::Id::Godunov))->getDeviceDataPtr());
        drKrnl.Qext = (entry.get(inner_keys::Wp::Id::FluxTarget))->getDeviceDataPtr();

        // the lift keeps a temporary where the face carries it per point
        real* tmpMem = nullptr;
        if constexpr (seissol::dynamicRupture::kernel::gpu_nodalFlux::TmpMaxMemRequiredInBytes >
                      0) {
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
    });
  }

  // the weight-adjusted cells have all of their fluxes in their buffers now, and add the weight
  // applied to them to their extended fields, which are cleared above (WeightAdjustable)
  weightAdjustment_.applyBatched(table, runtime, [](auto& krnl, auto& entry) {
    krnl.Qext = (entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr();
  });

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
