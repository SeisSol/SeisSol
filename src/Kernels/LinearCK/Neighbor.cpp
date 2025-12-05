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
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Parallel/Runtime/Stream.h"

#include <Common/Constants.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/Typedefs.h>
#include <Memory/Descriptor/LTS.h>
#include <Parallel/Runtime/Stream.h>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <stdint.h>

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#endif

#include <utils/logger.h>

#ifndef NDEBUG
#include "Alignment.h"
#endif

#include "Memory/GlobalData.h"

namespace seissol::kernels::solver::linearck {

template <typename Cfg>
void Neighbor<Cfg>::setGlobalData(const GlobalData& global) {

  m_nfKrnlPrototype.rDivM = global.get<Cfg>().changeOfBasisMatrices;
  m_nfKrnlPrototype.rT = global.get<Cfg>().neighborChangeOfBasisMatricesTransposed;
  m_nfKrnlPrototype.fP = global.get<Cfg>().neighborFluxMatrices;
  m_drKrnlPrototype.V3mTo2nTWDivM = global.get<Cfg>().nodalFluxMatrices;

#ifdef ACL_DEVICE
  assert(global.onDevice != nullptr);

#ifdef USE_PREMULTIPLY_FLUX
  deviceNfKrnlPrototype.minusFluxMatrices = global.get<Cfg, Executor::Device>().minusFluxMatrices;
#else
  deviceNfKrnlPrototype.rDivM = global.get<Cfg, Executor::Device>().changeOfBasisMatrices;
  deviceNfKrnlPrototype.rT =
      global.get<Cfg, Executor::Device>().neighborChangeOfBasisMatricesTransposed;
  deviceNfKrnlPrototype.fP = global.get<Cfg, Executor::Device>().neighborFluxMatrices;
#endif
  deviceDrKrnlPrototype.V3mTo2nTWDivM = global.get<Cfg, Executor::Device>().nodalFluxMatrices;
#endif
}

template <typename Cfg>
void Neighbor<Cfg>::computeNeighborsIntegral(LTS::Ref<Cfg>& data,
                                             const CellDRMapping<Cfg> (&cellDrMapping)[4],
                                             real* timeIntegrated[4],
                                             real* faceNeighborsPrefetch[4]) {
  assert(reinterpret_cast<uintptr_t>(data.template get<LTS::Dofs>()) % Alignment == 0);

  for (std::size_t face = 0; face < Cell::NumFaces; face++) {
    switch (data.template get<LTS::CellInformation>().faceTypes[face]) {
    case FaceType::Regular:
      // Fallthrough intended
    case FaceType::Periodic: {
      // Standard neighboring flux
      // Compute the neighboring elements flux matrix id.
      assert(reinterpret_cast<uintptr_t>(timeIntegrated[face]) % Alignment == 0);
      assert(data.template get<LTS::CellInformation>().faceRelations[face][0] < Cell::NumFaces &&
             data.template get<LTS::CellInformation>().faceRelations[face][1] < 3);
      kernel::neighboringFlux<Cfg> nfKrnl = m_nfKrnlPrototype;
      nfKrnl.Q = data.template get<LTS::Dofs>();
      nfKrnl.I = timeIntegrated[face];
      nfKrnl.AminusT = data.template get<LTS::NeighboringIntegration>().nAmNm1[face];
      nfKrnl._prefetch.I = faceNeighborsPrefetch[face];
      nfKrnl.execute(data.template get<LTS::CellInformation>().faceRelations[face][1],
                     data.template get<LTS::CellInformation>().faceRelations[face][0],
                     face);
      break;
    }
    case FaceType::DynamicRupture: {
      // No neighboring cell contribution, interior bc.
      assert(reinterpret_cast<uintptr_t>(cellDrMapping[face].godunov) % Alignment == 0);

      dynamicRupture::kernel::nodalFlux<Cfg> drKrnl = m_drKrnlPrototype;
      drKrnl.fluxSolver = cellDrMapping[face].fluxSolver;
      drKrnl.QInterpolated = cellDrMapping[face].godunov;
      drKrnl.Q = data.template get<LTS::Dofs>();
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

template <typename Cfg>
void Neighbor<Cfg>::computeBatchedNeighborsIntegral(
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& table,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  using namespace seissol::recording;
  kernel::gpu_neighboringFlux<Cfg> neighFluxKrnl = deviceNfKrnlPrototype;
  dynamicRupture::kernel::gpu_nodalFlux<Cfg> drKrnl = deviceDrKrnlPrototype;

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

template <typename Cfg>
void Neighbor<Cfg>::flopsNeighborsIntegral(
    const std::array<FaceType, Cell::NumFaces>& faceTypes,
    const std::array<std::array<uint8_t, 2>, Cell::NumFaces>& neighboringIndices,
    const CellDRMapping<Cfg> (&cellDrMapping)[4],
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
      nonZeroFlops += kernel::neighboringFlux<Cfg>::nonZeroFlops(
          neighboringIndices[face][1], neighboringIndices[face][0], face);
      hardwareFlops += kernel::neighboringFlux<Cfg>::hardwareFlops(
          neighboringIndices[face][1], neighboringIndices[face][0], face);
      break;
    case FaceType::DynamicRupture:
      drNonZeroFlops += dynamicRupture::kernel::nodalFlux<Cfg>::nonZeroFlops(
          cellDrMapping[face].side, cellDrMapping[face].faceRelation);
      drHardwareFlops += dynamicRupture::kernel::nodalFlux<Cfg>::hardwareFlops(
          cellDrMapping[face].side, cellDrMapping[face].faceRelation);
      break;
    default:
      // Handled in local kernel
      break;
    }
  }
}

template <typename Cfg>
std::uint64_t Neighbor<Cfg>::bytesNeighborsIntegral() {
  std::uint64_t reals = 0;

  // 4 * tElasticDOFS load, DOFs load, DOFs write
  reals += 4 * tensor::I<Cfg>::size() + 2 * tensor::Q<Cfg>::size();
  // flux solvers load
  reals += static_cast<std::uint64_t>(4 * tensor::AminusT<Cfg>::size());

  return reals * sizeof(real);
}

#define SEISSOL_CONFIGITER(cfg) template class Neighbor<cfg>;
#include "ConfigIncludeLinearCK.h"

#define SEISSOL_CONFIGITER(cfg) template class Neighbor<cfg>;
#include "ConfigIncludeSTP.h"

} // namespace seissol::kernels::solver::linearck
