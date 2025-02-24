// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
// SPDX-FileCopyrightText: 2014-2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#include "Kernels/Neighbor.h"

#include <Common/Constants.h>
#include <DataTypes/ConditionalTable.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Interface.h>
#include <Kernels/NeighborBase.h>
#include <Kernels/Precision.h>
#include <Parallel/Runtime/Stream.h>
#include <cassert>
#include <cstddef>
#include <generated_code/tensor.h>
#include <stdint.h>

#include "utils/logger.h"

namespace seissol::kernels {

void NeighborBase::checkGlobalData(const GlobalData* global, size_t alignment) {
#ifndef NDEBUG
  for (int neighbor = 0; neighbor < 4; ++neighbor) {
    assert((reinterpret_cast<uintptr_t>(global->changeOfBasisMatrices(neighbor))) % alignment == 0);
    assert((reinterpret_cast<uintptr_t>(global->localChangeOfBasisMatricesTransposed(neighbor))) %
               alignment ==
           0);
    assert(
        (reinterpret_cast<uintptr_t>(global->neighbourChangeOfBasisMatricesTransposed(neighbor))) %
            alignment ==
        0);
  }

  for (int h = 0; h < 3; ++h) {
    assert((reinterpret_cast<uintptr_t>(global->neighbourFluxMatrices(h))) % alignment == 0);
  }

  for (int i = 0; i < 4; ++i) {
    for (int h = 0; h < 3; ++h) {
      assert((reinterpret_cast<uintptr_t>(global->nodalFluxMatrices(i, h))) % alignment == 0);
    }
  }
#endif
}

void Neighbor::setHostGlobalData(const GlobalData* global) {
  checkGlobalData(global, Alignment);
  m_nfKrnlPrototype.rDivM = global->changeOfBasisMatrices;
  m_nfKrnlPrototype.rT = global->neighbourChangeOfBasisMatricesTransposed;
  m_nfKrnlPrototype.fP = global->neighbourFluxMatrices;
  m_drKrnlPrototype.V3mTo2nTWDivM = global->nodalFluxMatrices;
}

void Neighbor::setGlobalData(const CompoundGlobalData& global) {
  setHostGlobalData(global.onHost);

#ifdef ACL_DEVICE
  assert(global.onDevice != nullptr);
  const auto deviceAlignment = device.api->getGlobMemAlignment();
  checkGlobalData(global.onDevice, deviceAlignment);

#ifdef USE_PREMULTIPLY_FLUX
  deviceNfKrnlPrototype.minusFluxMatrices = global.onDevice->minusFluxMatrices;
#else
  deviceNfKrnlPrototype.rDivM = global.onDevice->changeOfBasisMatrices;
  deviceNfKrnlPrototype.rT = global.onDevice->neighbourChangeOfBasisMatricesTransposed;
  deviceNfKrnlPrototype.fP = global.onDevice->neighbourFluxMatrices;
#endif
  deviceDrKrnlPrototype.V3mTo2nTWDivM = global.onDevice->nodalFluxMatrices;
#endif
}

void Neighbor::computeNeighborsIntegral(NeighborData& data,
                                        const CellDRMapping (&cellDrMapping)[4],
                                        real* timeIntegrated[4],
                                        real* faceNeighborsPrefetch[4]) {
  assert(reinterpret_cast<uintptr_t>(data.dofs()) % Alignment == 0);

  for (unsigned int face = 0; face < 4; face++) {
    switch (data.cellInformation().faceTypes[face]) {
    case FaceType::Regular:
      // Fallthrough intended
    case FaceType::Periodic: {
      // Standard neighboring flux
      // Compute the neighboring elements flux matrix id.
      assert(reinterpret_cast<uintptr_t>(timeIntegrated[face]) % Alignment == 0);
      assert(data.cellInformation().faceRelations[face][0] < 4 &&
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

  real* tmpMem = nullptr;
  auto resetDeviceCurrentState = [this](size_t counter) {
    for (size_t i = 0; i < counter; ++i) {
      this->device.api->popStackMemory();
    }
  };

  for (size_t face = 0; face < 4; face++) {
    size_t streamCounter{0};

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
                  (entry.get(inner_keys::Wp::Id::AminusT))->getDeviceDataPtr());

              tmpMem = reinterpret_cast<real*>(
                  device.api->getStackMemory(neighFluxKrnl.TmpMaxMemRequiredInBytes * numElements));
              neighFluxKrnl.linearAllocator.initialize(tmpMem);

              neighFluxKrnl.streamPtr = stream;
              (neighFluxKrnl.*neighFluxKrnl.ExecutePtrs[faceRelation])();
              ++streamCounter;
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

              tmpMem = reinterpret_cast<real*>(
                  device.api->getStackMemory(drKrnl.TmpMaxMemRequiredInBytes * numElements));
              drKrnl.linearAllocator.initialize(tmpMem);

              drKrnl.streamPtr = stream;
              (drKrnl.*drKrnl.ExecutePtrs[faceRelation])();
              ++streamCounter;
            }
          }
        });

    resetDeviceCurrentState(streamCounter);
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

void Neighbor::flopsNeighborsIntegral(const FaceType faceTypes[4],
                                      const int neighboringIndices[4][2],
                                      const CellDRMapping (&cellDrMapping)[4],
                                      unsigned int& nonZeroFlops,
                                      unsigned int& hardwareFlops,
                                      long long& drNonZeroFlops,
                                      long long& drHardwareFlops) {
  // reset flops
  nonZeroFlops = 0;
  hardwareFlops = 0;
  drNonZeroFlops = 0;
  drHardwareFlops = 0;

  for (unsigned int face = 0; face < 4; face++) {
    // compute the neighboring elements flux matrix id.
    switch (faceTypes[face]) {
    case FaceType::Regular:
      // Fallthrough intended
    case FaceType::Periodic:
      // regular neighbor
      assert(neighboringIndices[face][0] < 4 && neighboringIndices[face][1] < 3);
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

unsigned Neighbor::bytesNeighborsIntegral() {
  unsigned reals = 0;

  // 4 * tElasticDOFS load, DOFs load, DOFs write
  reals += 4 * tensor::I::size() + 2 * tensor::Q::size();
  // flux solvers load
  reals += 4 * tensor::AminusT::size();

  return reals * sizeof(real);
}

} // namespace seissol::kernels
