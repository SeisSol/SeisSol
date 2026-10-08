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
#include "Config.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/LTS.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <stdint.h>
#include <utility>

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#endif

#ifndef ACL_DEVICE
#include <utils/logger.h>
#endif

#ifndef NDEBUG
#include "Alignment.h"
#endif

namespace seissol::kernels::solver::linearck {

template <typename Cfg>
void Neighbor<Cfg>::setGlobalData(const CompoundGlobalData<Cfg>& global) {

  nfKrnlPrototype_.bindGlobals(*global.onHost);
  drKrnlPrototype_.bindGlobals(*global.onHost);

#ifdef ACL_DEVICE
  assert(global.onDevice != nullptr);

  deviceNfKrnlPrototype_.bindGlobals(*global.onDevice);
  deviceDrKrnlPrototype_.bindGlobals(*global.onDevice);
#endif
}

template <typename Cfg>
void Neighbor<Cfg>::computeNeighborsIntegral(
    LTS::Ref<Cfg>& data,
    const std::array<real*, Cell::NumFaces>& timeIntegrated,
    const std::array<real*, Cell::NumFaces>& faceNeighborsPrefetch) {
  // The neighbouring flux family is indexed by the neighbouring side and the own face. The face
  // orientation index is not part of it, since the canonical vertex numbering pins it to zero on
  // every interior face.
  static_assert(std::size(kernel::neighboringFlux<Cfg>::ExecutePtrs) ==
                Cell::NumFaces * Cell::NumFaces);

  assert(reinterpret_cast<uintptr_t>(data.template get<LTS::Dofs>()) % Vectorsize == 0);
  const auto& cellDrMapping = data.template get<LTS::DRMapping>();

  for (std::size_t face = 0; face < Cell::NumFaces; face++) {
    switch (data.template get<LTS::CellInformation>().faceTypes[face]) {
    case FaceType::Regular: {
      // Standard neighboring flux
      // Compute the neighboring elements flux matrix id.
      assert(reinterpret_cast<uintptr_t>(timeIntegrated[face]) % Vectorsize == 0);
      assert(data.template get<LTS::CellInformation>().faceRelations[face][0] < Cell::NumFaces &&
             data.template get<LTS::CellInformation>().faceRelations[face][1] == 0);
      kernel::neighboringFlux<Cfg> nfKrnl = nfKrnlPrototype_;
      nfKrnl.Q = data.template get<LTS::Dofs>();
      nfKrnl.I = timeIntegrated[face];
      nfKrnl.AminusT = data.template get<LTS::NeighboringIntegration>().nAmNm1[face];
      nfKrnl._prefetch.I = faceNeighborsPrefetch[face];
      nfKrnl.execute(data.template get<LTS::CellInformation>().faceRelations[face][0], face);
      break;
    }
    case FaceType::DynamicRupture: {
      // No neighboring cell contribution, interior bc.
      assert(reinterpret_cast<uintptr_t>(cellDrMapping[face].godunov) % Vectorsize == 0);

      dynamicRupture::kernel::nodalFlux<Cfg> drKrnl = drKrnlPrototype_;
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
  static_assert(std::size(kernel::gpu_neighboringFlux<Cfg>::ExecutePtrs) ==
                *seissol::recording::FaceRelations::Count);
  static_assert(std::size(dynamicRupture::kernel::gpu_nodalFlux<Cfg>::ExecutePtrs) ==
                *seissol::recording::DrFaceRelations::Count);

  using namespace seissol::recording;
  kernel::gpu_neighboringFlux<Cfg> neighFluxKrnl = deviceNfKrnlPrototype_;
  dynamicRupture::kernel::gpu_nodalFlux<Cfg> drKrnl = deviceDrKrnlPrototype_;

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    runtime.envMany(
        (*FaceRelations::PerFace) + (*DrFaceRelations::PerFace), [&](void* stream, size_t i) {
          if (i < (*FaceRelations::PerFace)) {
            // regular and periodic
            const auto faceRelation = i + (*FaceRelations::PerFace) * face;

            const ConditionalKey key(
                *KernelNames::NeighborFlux, *FaceKinds::Regular, face, faceRelation);

            if (table.find(key) != table.end()) {
              auto& entry = table[key];

              const auto numElements = (entry.get<real*>(inner_keys::Wp::Id::Dofs))->getSize();
              neighFluxKrnl.numElements = numElements;

              neighFluxKrnl.Q = (entry.get<real*>(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
              neighFluxKrnl.I = const_cast<const real**>(
                  (entry.get<real*>(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());
              neighFluxKrnl.AminusT = const_cast<const real**>(
                  entry.get<real*>(inner_keys::Wp::Id::NeighborIntegrationData)
                      ->getDeviceDataPtr());

              SEISSOL_ARRAY_OFFSET_ASSERT(NeighboringIntegrationData<Cfg>, nAmNm1);
              neighFluxKrnl.extraOffset_AminusT =
                  SEISSOL_ARRAY_OFFSET(NeighboringIntegrationData<Cfg>, nAmNm1, face);

              real* tmpMem = reinterpret_cast<real*>(device_.api().allocMemAsync(
                  seissol::kernel::gpu_neighboringFlux<Cfg>::TmpMaxMemRequiredInBytes * numElements,
                  stream));
              neighFluxKrnl.linearAllocator.initialize(tmpMem);

              neighFluxKrnl.streamPtr = stream;
              (neighFluxKrnl.*
               seissol::kernel::gpu_neighboringFlux<Cfg>::ExecutePtrs[faceRelation])();
              device_.api().freeMemAsync(reinterpret_cast<void*>(tmpMem), stream);
            }
          } else {
            // the side is the minor index here, cf. the NeighIntegrationRecorder
            const auto faceRelation = face + Cell::NumFaces * (i - (*FaceRelations::PerFace));

            const ConditionalKey key(
                *KernelNames::NeighborFlux, *FaceKinds::DynamicRupture, face, faceRelation);

            if (table.find(key) != table.end()) {
              auto& entry = table[key];

              const auto numElements = (entry.get<real*>(inner_keys::Wp::Id::Dofs))->getSize();
              drKrnl.numElements = numElements;

              drKrnl.fluxSolver = const_cast<const real**>(
                  (entry.get<real*>(inner_keys::Wp::Id::FluxSolver))->getDeviceDataPtr());
              drKrnl.QInterpolated = const_cast<const real**>(
                  (entry.get<real*>(inner_keys::Wp::Id::Godunov))->getDeviceDataPtr());
              drKrnl.Q = (entry.get<real*>(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();

              real* tmpMem = reinterpret_cast<real*>(device_.api().allocMemAsync(
                  seissol::dynamicRupture::kernel::gpu_nodalFlux<Cfg>::TmpMaxMemRequiredInBytes *
                      numElements,
                  stream));
              drKrnl.linearAllocator.initialize(tmpMem);

              drKrnl.streamPtr = stream;
              (drKrnl.*
               seissol::dynamicRupture::kernel::gpu_nodalFlux<Cfg>::ExecutePtrs[faceRelation])();
              device_.api().freeMemAsync(reinterpret_cast<void*>(tmpMem), stream);
            }
          }
        });
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

template <typename Cfg>
std::pair<PerformanceEstimate, PerformanceEstimate> Neighbor<Cfg>::metrics(
    const std::array<FaceType, Cell::NumFaces>& faceTypes,
    const std::array<std::array<uint8_t, 2>, Cell::NumFaces>& neighboringIndices,
    const std::array<CellDRMapping<Cfg>, Cell::NumFaces>& cellDrMapping) const {
  // reset flops
  PerformanceEstimate neigh;
  PerformanceEstimate neighDR;

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // compute the neighboring elements flux matrix id.
    switch (faceTypes[face]) {
    case FaceType::Regular:
      // regular neighbor
      assert(neighboringIndices[face][0] < Cell::NumFaces && neighboringIndices[face][1] == 0);
      neigh += PerformanceEstimate::fromKernel<kernel::neighboringFlux<Cfg>>(
          neighboringIndices[face][0], face);
      break;
    case FaceType::DynamicRupture:
      neighDR += PerformanceEstimate::fromKernel<dynamicRupture::kernel::nodalFlux<Cfg>>(
          cellDrMapping[face].side, cellDrMapping[face].faceRelation);
      break;
    default:
      // Handled in local kernel
      break;
    }
  }

  // legacy memory estimate
  std::uint64_t reals = 0;

  // 4 * tElasticDOFS load, DOFs load, DOFs write
  reals += 4 * tensor::I<Cfg>::size() + 2 * tensor::Q<Cfg>::size();
  // flux solvers load
  reals += static_cast<std::uint64_t>(4 * tensor::AminusT<Cfg>::size());

  neigh.bytes = reals * sizeof(real);

  return {neigh, neighDR};
}

#define SEISSOL_INSTANTIATE(Cfg) template class Neighbor<Cfg>;
SEISSOL_FOR_EACH_CONFIG_LINEARCK(SEISSOL_INSTANTIATE)
SEISSOL_FOR_EACH_CONFIG_STP(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::kernels::solver::linearck
