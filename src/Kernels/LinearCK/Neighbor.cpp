// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2014-2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#include "Kernels/LinearCK/Neighbor.h"

#include "Common/Constants.h"
#include "Common/Marker.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/Typedefs.h"
#include "Kernels/FaceTraces.h"
#include "Kernels/Precision.h"
#include "Kernels/StarOperands.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Model/OperatorLayout.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <stdint.h>
#include <utility>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#endif

#ifndef NDEBUG
#include "Alignment.h"
#endif

namespace seissol::kernels::solver::linearck {

// The corrector family is indexed by the face and by the side of it the cell is on, which says in
// whose parametrisation the traces of the face are written. The face orientation index is not part
// of it, since the canonical vertex numbering pins it to zero on every interior face.
static_assert(std::size(kernel::faceCorrector::ExecutePtrs) == Cell::NumFaces * 2);

#ifdef ACL_DEVICE
static_assert(std::size(kernel::gpu_faceCorrector::ExecutePtrs) == Cell::NumFaces * 2);
static_assert(std::size(dynamicRupture::kernel::gpu_nodalFlux::ExecutePtrs) ==
              *seissol::recording::DrFaceRelations::Count);
#endif

void Neighbor::setGlobalData(const CompoundGlobalData& global) {

  correctorKrnlPrototype_.bindGlobals(*global.onHost);
  drKrnlPrototype_.bindGlobals(*global.onHost);
  weightAdjustment_.setGlobalData(global);

#ifdef ACL_DEVICE
  assert(global.onDevice != nullptr);

  deviceCorrectorKrnlPrototype_.bindGlobals(*global.onDevice);
  deviceDrKrnlPrototype_.bindGlobals(*global.onDevice);
#endif
}

void Neighbor::computeNeighborsIntegral(
    LTS::Ref& data, const std::array<const real*, Cell::NumFaces>& neighborTraces) {
  assert(reinterpret_cast<uintptr_t>(data.get<LTS::Dofs>()) % Vectorsize == 0);
  const auto& cellDrMapping = data.get<LTS::DRMapping>();

  // a weight-adjusted cell collects what its faces lift in its flux buffer, which the local
  // integral started, and applies the inverse weight once they are all in (WeightAdjustable)
  const bool weightAdjusted = WeightAdjustable && data.get<LTS::CellInformation>().weightAdjusted;
  real* fluxTarget = weightAdjusted ? data.get<LTS::FluxBuffer>() : data.get<LTS::Dofs>();

  kernel::faceCorrector correctorKrnl = correctorKrnlPrototype_;
  correctorKrnl.Q = fluxTarget;
  kernels::applyFaceCorrectors(correctorKrnl, neighborTraces, data);

  // the flux the fault supplies on its faces
  for (std::size_t face = 0; face < Cell::NumFaces; face++) {
    if (data.get<LTS::CellInformation>().faceTypes[face] == FaceType::DynamicRupture) {
      assert(reinterpret_cast<uintptr_t>(cellDrMapping[face].godunov) % Vectorsize == 0);

      dynamicRupture::kernel::nodalFlux drKrnl = drKrnlPrototype_;
      kernels::bindFaultFluxOperands(drKrnl, cellDrMapping[face].fluxSolver);
      drKrnl.QInterpolated = cellDrMapping[face].godunov;
      drKrnl.Q = fluxTarget;
      drKrnl.execute(cellDrMapping[face].side, cellDrMapping[face].faceRelation);
    }
  }

  if (weightAdjusted) {
    weightAdjustment_.apply(data, [&](auto& krnl) { krnl.Q = data.get<LTS::Dofs>(); });
  }
}

void Neighbor::computeBatchedNeighborsIntegral(
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& table,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  using namespace seissol::recording;

  // the corrector of every face, into what the faces of the cells lift into (FluxTarget)
  kernels::applyFaceCorrectorsBatched(
      deviceCorrectorKrnlPrototype_, table, runtime, [](auto& krnl, auto& entry) {
        krnl.Q = (entry.get(inner_keys::Wp::Id::FluxTarget))->getDeviceDataPtr();
      });

  // the flux the fault supplies on its faces; the face relations of one face are different cells
  dynamicRupture::kernel::gpu_nodalFlux drKrnl = deviceDrKrnlPrototype_;
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    runtime.envMany((*DrFaceRelations::Count), [&](void* stream, size_t faceRelation) {
      const ConditionalKey key(
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
        drKrnl.Q = (entry.get(inner_keys::Wp::Id::FluxTarget))->getDeviceDataPtr();

        real* tmpMem = reinterpret_cast<real*>(device_.api->allocMemAsync(
            seissol::dynamicRupture::kernel::gpu_nodalFlux::TmpMaxMemRequiredInBytes * numElements,
            stream));
        drKrnl.linearAllocator.initialize(tmpMem);

        drKrnl.streamPtr = stream;
        (drKrnl.*seissol::dynamicRupture::kernel::gpu_nodalFlux::ExecutePtrs[faceRelation])();
        device_.api->freeMemAsync(reinterpret_cast<void*>(tmpMem), stream);
      }
    });
  }

  // the weight-adjusted cells have all of their fluxes in their buffers now (WeightAdjustable)
  weightAdjustment_.applyBatched(table, runtime, [](auto& krnl, auto& entry) {
    krnl.Q = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
  });
#else
  logError() << "No GPU implementation provided";
#endif
}

std::pair<PerformanceEstimate, PerformanceEstimate>
    Neighbor::metrics(const std::array<FaceType, Cell::NumFaces>& faceTypes,
                      const std::array<CellDRMapping, Cell::NumFaces>& cellDrMapping) const {
  // reset flops
  PerformanceEstimate neigh;
  PerformanceEstimate neighDR;

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // the corrector of the face; both sides of a face cost the same
    if (appliesFaceCorrector(faceTypes[face])) {
      neigh += PerformanceEstimate::fromKernel<kernel::faceCorrector>(face, 0);
    }
    if (faceTypes[face] == FaceType::DynamicRupture) {
      neighDR += PerformanceEstimate::fromKernel<dynamicRupture::kernel::nodalFlux>(
          cellDrMapping[face].side, cellDrMapping[face].faceRelation);
    }
  }

  // legacy memory estimate
  std::uint64_t reals = 0;

  // 8 traces load, DOFs load, DOFs write
  reals += 2 * Cell::NumFaces * tensor::faceTrace::size() + 2 * tensor::Q::size();

  neigh.bytes = reals * sizeof(real);

  return {neigh, neighDR};
}

} // namespace seissol::kernels::solver::linearck
