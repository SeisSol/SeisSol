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

#include "Kernels/LinearCK/NeighborBase.h"

#include <Common/Constants.h>
#include <DataTypes/ConditionalTable.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Interface.h>
#include <Kernels/Precision.h>
#include <Parallel/Runtime/Stream.h>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <generated_code/tensor.h>
#include <stdint.h>

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#endif

#include "utils/logger.h"

#ifndef NDEBUG
#include "Alignment.h"
#endif

namespace seissol::kernels::solver::linearck {
void Neighbor::setGlobalData(const CompoundGlobalData& global) {

  m_nfKrnlPrototype.rDivM = global.onHost->changeOfBasisMatrices;
  m_nfKrnlPrototype.rT = global.onHost->neighborChangeOfBasisMatricesTransposed;
  m_nfKrnlPrototype.fP = global.onHost->neighborFluxMatrices;
  m_drKrnlPrototype.V3mTo2nTWDivM = global.onHost->nodalFluxMatrices;

#ifdef ACL_DEVICE
  assert(global.onDevice != nullptr);
  const auto deviceAlignment = device.api->getGlobMemAlignment();

#ifdef USE_PREMULTIPLY_FLUX
  deviceNfKrnlPrototype.minusFluxMatrices = global.onDevice->minusFluxMatrices;
#else
  deviceNfKrnlPrototype.rDivM = global.onDevice->changeOfBasisMatrices;
  deviceNfKrnlPrototype.rT = global.onDevice->neighborChangeOfBasisMatricesTransposed;
  deviceNfKrnlPrototype.fP = global.onDevice->neighborFluxMatrices;
#endif
  deviceDrKrnlPrototype.V3mTo2nTWDivM = global.onDevice->nodalFluxMatrices;
#endif
}

void Neighbor::computeNeighborsIntegral(NeighborData& data,
                                        const CellDRMapping (&cellDrMapping)[4],
                                        real* timeIntegrated[4],
                                        real* faceNeighborsPrefetch[4]) {
  assert(reinterpret_cast<uintptr_t>(data.dofs()) % Alignment == 0);

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    switch (data.cellInformation().faceTypes[face]) {
    case FaceType::Regular:
      // Fallthrough intended
    case FaceType::Periodic: {
      // Standard neighboring flux
      // Compute the neighboring elements flux matrix id.
      assert(reinterpret_cast<uintptr_t>(timeIntegrated[face]) % Alignment == 0);
      assert(data.cellInformation().faceRelations[face][0] < Cell::NumFaces &&
             data.cellInformation().faceRelations[face][1] < 3);
      kernel::neighboringFlux nfKrnl = m_nfKrnlPrototype;
      nfKrnl.Q = data.dofs();
      nfKrnl.I = timeIntegrated[face];
      nfKrnl.AminusT = data.neighboringIntegration().nAmNm1[face];
      nfKrnl._prefetch.I = faceNeighborsPrefetch[face];
      nfKrnl.execute(data.cellInformation().faceRelations[face][1],
                     data.cellInformation().faceRelations[face][0],
                     face);
      break;
    }
    case FaceType::DynamicRupture: {
      // No neighboring cell contribution, interior bc.
      assert(reinterpret_cast<uintptr_t>(cellDrMapping[face].godunov) % Alignment == 0);

      dynamicRupture::kernel::nodalFlux drKrnl = m_drKrnlPrototype;
      drKrnl.fluxSolver = cellDrMapping[face].fluxSolver;
      drKrnl.QInterpolated = cellDrMapping[face].godunov;
      drKrnl.Q = data.dofs();
      drKrnl._prefetch.I = faceNeighborsPrefetch[face];
      drKrnl.execute(cellDrMapping[face].side, cellDrMapping[face].faceRelation);
      break;
    }
    default:
      // No contribution for all other cases.
      // Note: some other bcs are handled in the local kernel.
      break;
    }
  }
}

void Neighbor::computeBatchedNeighborsIntegral(ConditionalPointersToRealsTable& table,
                                               seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  kernel::gpu_neighboringFlux neighFluxKrnl = deviceNfKrnlPrototype;
  dynamicRupture::kernel::gpu_nodalFlux drKrnl = deviceDrKrnlPrototype;

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    runtime.envMany(
        (*FaceRelations::Count) + (*DrFaceRelations::Count), [&](void* stream, size_t i) {
          if (i < (*FaceRelations::Count)) {
            // regular and periodic
            unsigned faceRelation = i;

            ConditionalKey key(*KernelNames::NeighborFlux,
                               (FaceKinds::Regular || FaceKinds::Periodic),
                               face,
                               faceRelation);

            if (table.find(key) != table.end()) {
              auto& entry = table[key];

              const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
              neighFluxKrnl.numElements = numElements;

              neighFluxKrnl.Q = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
              neighFluxKrnl.I = const_cast<const real**>(
                  (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());
              neighFluxKrnl.AminusT = const_cast<const real**>(
                  entry.get(inner_keys::Wp::Id::NeighborIntegrationData)->getDeviceDataPtr());
              neighFluxKrnl.extraOffset_AminusT =
                  SEISSOL_ARRAY_OFFSET(NeighboringIntegrationData, nAmNm1, face);

              real* tmpMem = reinterpret_cast<real*>(device.api->allocMemAsync(
                  neighFluxKrnl.TmpMaxMemRequiredInBytes * numElements, stream));
              neighFluxKrnl.linearAllocator.initialize(tmpMem);

              neighFluxKrnl.streamPtr = stream;
              (neighFluxKrnl.*neighFluxKrnl.ExecutePtrs[faceRelation])();
              device.api->freeMemAsync(reinterpret_cast<void*>(tmpMem), stream);
            }
          } else {
            unsigned faceRelation = i - (*FaceRelations::Count);

            ConditionalKey key(
                *KernelNames::NeighborFlux, *FaceKinds::DynamicRupture, face, faceRelation);

            if (table.find(key) != table.end()) {
              auto& entry = table[key];

              const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
              drKrnl.numElements = numElements;

              drKrnl.fluxSolver = const_cast<const real**>(
                  (entry.get(inner_keys::Wp::Id::FluxSolver))->getDeviceDataPtr());
              drKrnl.QInterpolated = const_cast<const real**>(
                  (entry.get(inner_keys::Wp::Id::Godunov))->getDeviceDataPtr());
              drKrnl.Q = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();

              real* tmpMem = reinterpret_cast<real*>(
                  device.api->allocMemAsync(drKrnl.TmpMaxMemRequiredInBytes * numElements, stream));
              drKrnl.linearAllocator.initialize(tmpMem);

              drKrnl.streamPtr = stream;
              (drKrnl.*drKrnl.ExecutePtrs[faceRelation])();
              device.api->freeMemAsync(reinterpret_cast<void*>(tmpMem), stream);
            }
          }
        });
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

void Neighbor::flopsNeighborsIntegral(const FaceType faceTypes[4],
                                      const int neighboringIndices[4][2],
                                      const CellDRMapping (&cellDrMapping)[4],
                                      std::uint64_t& nonZeroFlops,
                                      std::uint64_t& hardwareFlops,
                                      std::uint64_t& drNonZeroFlops,
                                      std::uint64_t& drHardwareFlops) {
  // reset flops
  nonZeroFlops = 0;
  hardwareFlops = 0;
  drNonZeroFlops = 0;
  drHardwareFlops = 0;

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // compute the neighboring elements flux matrix id.
    switch (faceTypes[face]) {
    case FaceType::Regular:
      // Fallthrough intended
    case FaceType::Periodic:
      // regular neighbor
      assert(neighboringIndices[face][0] < Cell::NumFaces && neighboringIndices[face][1] < 3);
      nonZeroFlops += kernel::neighboringFlux::nonZeroFlops(
          neighboringIndices[face][1], neighboringIndices[face][0], face);
      hardwareFlops += kernel::neighboringFlux::hardwareFlops(
          neighboringIndices[face][1], neighboringIndices[face][0], face);
      break;
    case FaceType::DynamicRupture:
      drNonZeroFlops += dynamicRupture::kernel::nodalFlux::nonZeroFlops(
          cellDrMapping[face].side, cellDrMapping[face].faceRelation);
      drHardwareFlops += dynamicRupture::kernel::nodalFlux::hardwareFlops(
          cellDrMapping[face].side, cellDrMapping[face].faceRelation);
      break;
    default:
      // Handled in local kernel
      break;
    }
  }
}

std::uint64_t Neighbor::bytesNeighborsIntegral() {
  std::uint64_t reals = 0;

  // 4 * tElasticDOFS load, DOFs load, DOFs write
  reals += 4 * tensor::I::size() + 2 * tensor::Q::size();
  // flux solvers load
  reals += static_cast<std::uint64_t>(4 * tensor::AminusT::size());

  return reals * sizeof(real);
}

} // namespace seissol::kernels::solver::linearck
