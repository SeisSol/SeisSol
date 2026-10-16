// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#include "Neighbor.h"

#include "Alignment.h"
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
#include <cstring>
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

namespace seissol::kernels::solver::linearckanelastic {

// The neighbouring flux family is indexed by the neighbouring side and the own face. The face
// orientation index is not part of it, since the canonical vertex numbering pins it to zero on
// every interior face.
#define SEISSOL_CHECK_FLUX_FAMILY(Cfg)                                                             \
  static_assert(std::size(seissol::kernel::neighborFluxExt<Cfg>::ExecutePtrs) ==                   \
                Cell::NumFaces * Cell::NumFaces);
SEISSOL_FOR_EACH_CONFIG_LINEARCKANELASTIC(SEISSOL_CHECK_FLUX_FAMILY)
#undef SEISSOL_CHECK_FLUX_FAMILY

template <typename Cfg>
void Neighbor<Cfg>::setGlobalData(const CompoundGlobalData<Cfg>& global) {
  nfKrnlPrototype_.bindGlobals(*global.onHost);
  drKrnlPrototype_.bindGlobals(*global.onHost);

#ifdef ACL_DEVICE
  deviceNfKrnlPrototype_.bindGlobals(*global.onDevice);
  deviceDrKrnlPrototype_.bindGlobals(*global.onDevice);
#endif
}

template <typename Cfg>
void Neighbor<Cfg>::computeNeighborsIntegral(
    LTS::Ref<Cfg>& data,
    const std::array<real*, Cell::NumFaces>& timeIntegrated,
    const std::array<real*, Cell::NumFaces>& faceNeighborsPrefetch) {
#ifndef NDEBUG
  for (std::size_t neighbor = 0; neighbor < Cell::NumFaces; ++neighbor) {
    // alignment of the time integrated dofs (only for linear interior)
    if (data.template get<LTS::CellInformation>().faceTypes[neighbor] == FaceType::Regular) {
      assert((reinterpret_cast<uintptr_t>(timeIntegrated[neighbor])) % Vectorsize == 0);
    }
  }
#endif

  const auto& cellDrMapping = data.template get<LTS::DRMapping>();

  // alignment of the degrees of freedom
  assert((reinterpret_cast<uintptr_t>(data.template get<LTS::Dofs>())) % Vectorsize == 0);

  alignas(PagesizeStack) real qext[tensor::Qext<Cfg>::size()] = {};

  kernel::neighborFluxExt<Cfg> nfKrnl = nfKrnlPrototype_;
  nfKrnl.Qext = qext;

  // iterate over faces
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // neighboring cell contribution only for interior faces
    if (data.template get<LTS::CellInformation>().faceTypes[face] == FaceType::Regular) {
      assert(data.template get<LTS::CellInformation>().faceRelations[face][0] < Cell::NumFaces &&
             data.template get<LTS::CellInformation>().faceRelations[face][1] == 0);

      nfKrnl.I = timeIntegrated[face];
      nfKrnl.AminusT = data.template get<LTS::NeighboringIntegration>().nAmNm1[face];
      nfKrnl._prefetch.I = faceNeighborsPrefetch[face];
      nfKrnl.execute(data.template get<LTS::CellInformation>().faceRelations[face][0], face);
    } else if (data.template get<LTS::CellInformation>().faceTypes[face] ==
               FaceType::DynamicRupture) {
      assert((reinterpret_cast<uintptr_t>(cellDrMapping[face].godunov)) % Vectorsize == 0);

      dynamicRupture::kernel::nodalFlux<Cfg> drKrnl = drKrnlPrototype_;
      drKrnl.fluxSolver = cellDrMapping[face].fluxSolver;
      drKrnl.QInterpolated = cellDrMapping[face].godunov;
      drKrnl.Qext = qext;
      drKrnl._prefetch.I = faceNeighborsPrefetch[face];
      drKrnl.execute(cellDrMapping[face].side, cellDrMapping[face].faceRelation);
    }
  }

  kernel::neighbor<Cfg> nKrnl = nKrnlPrototype_;
  nKrnl.Qext = qext;
  nKrnl.Q = data.template get<LTS::Dofs>();
  nKrnl.Qane = data.template get<LTS::DofsAne>();
  nKrnl.w = data.template get<LTS::NeighboringIntegration>().specific.w;

  nKrnl.execute();
}

template <typename Cfg>
std::pair<PerformanceEstimate, PerformanceEstimate> Neighbor<Cfg>::metrics(
    const std::array<FaceType, Cell::NumFaces>& faceTypes,
    const std::array<std::array<uint8_t, 2>, Cell::NumFaces>& neighboringIndices,
    const std::array<CellDRMapping<Cfg>, Cell::NumFaces>& cellDrMapping) const {

  PerformanceEstimate regular;
  PerformanceEstimate dr;

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // neighboring cell contribution only for interior faces
    if (faceTypes[face] == FaceType::Regular) {
      assert(neighboringIndices[face][0] < Cell::NumFaces && neighboringIndices[face][1] == 0);

      regular += PerformanceEstimate::fromKernel<seissol::kernel::neighborFluxExt<Cfg>>(
          neighboringIndices[face][0], face);
    } else if (faceTypes[face] == FaceType::DynamicRupture) {
      dr += PerformanceEstimate::fromKernel<dynamicRupture::kernel::nodalFlux<Cfg>>(
          cellDrMapping[face].side, cellDrMapping[face].faceRelation);
    }
  }

  regular += PerformanceEstimate::fromKernel<kernel::neighbor<Cfg>>();

  // legacy memory estimate
  std::uint64_t reals = 0;

  // 4 * tElasticDOFS load, DOFs load, DOFs write
  reals += 4 * tensor::I<Cfg>::size() + 2 * tensor::Q<Cfg>::size() + 2 * tensor::Qane<Cfg>::size();
  // flux solvers load
  reals += 4 * tensor::AminusT<Cfg>::size() + tensor::w<Cfg>::size();

  regular.bytes = reals * sizeof(real);

  return {regular, dr};
}

template <typename Cfg>
void Neighbor<Cfg>::computeBatchedNeighborsIntegral(
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& table,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE

  using namespace seissol::recording;
  kernel::gpu_neighborFluxExt<Cfg> neighFluxKrnl = deviceNfKrnlPrototype_;
  dynamicRupture::kernel::gpu_nodalFlux<Cfg> drKrnl = deviceDrKrnlPrototype_;

  {
    const ConditionalKey key(KernelNames::Time || KernelNames::Volume);
    if (table.find(key) != table.end()) {
      auto& entry = table[key];
      device::DeviceInstance::instance().algorithms().setToValue(
          (entry.get<real*>(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr(),
          static_cast<real>(0.0),
          tensor::Qext<Cfg>::Size,
          (entry.get<real*>(inner_keys::Wp::Id::DofsExt))->getSize(),
          runtime.stream());
    }
  }

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    runtime.envMany(
        (*FaceRelations::PerFace) + (*DrFaceRelations::PerFace), [&](void* stream, size_t i) {
          // regular and periodic
          if (i < (*FaceRelations::PerFace)) {
            // regular and periodic
            const auto faceRelation = i + (*FaceRelations::PerFace) * face;

            const ConditionalKey key(
                *KernelNames::NeighborFlux, *FaceKinds::Regular, face, faceRelation);

            if (table.find(key) != table.end()) {
              auto& entry = table[key];

              const auto numElements = (entry.get<real*>(inner_keys::Wp::Id::Dofs))->getSize();
              neighFluxKrnl.numElements = numElements;

              neighFluxKrnl.Qext =
                  (entry.get<real*>(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr();
              neighFluxKrnl.I = const_cast<const real**>(
                  (entry.get<real*>(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());
              neighFluxKrnl.AminusT = const_cast<const real**>(
                  entry.get<real*>(inner_keys::Wp::Id::NeighborIntegrationData)
                      ->getDeviceDataPtr());

              SEISSOL_ARRAY_OFFSET_ASSERT(NeighboringIntegrationData<Cfg>, nAmNm1);
              neighFluxKrnl.extraOffset_AminusT =
                  SEISSOL_ARRAY_OFFSET(NeighboringIntegrationData<Cfg>, nAmNm1, face);

              neighFluxKrnl.streamPtr = stream;
              (neighFluxKrnl.*
               seissol::kernel::gpu_neighborFluxExt<Cfg>::ExecutePtrs[faceRelation])();
            }
          } else {
            // Dynamic Rupture
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
              drKrnl.Qext = (entry.get<real*>(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr();

              drKrnl.streamPtr = stream;
              (drKrnl.*
               seissol::dynamicRupture::kernel::gpu_nodalFlux<Cfg>::ExecutePtrs[faceRelation])();
            }
          }
        });
  }

  const ConditionalKey key(KernelNames::Time || KernelNames::Volume);
  if (table.find(key) != table.end()) {
    auto& entry = table[key];
    kernel::gpu_neighbor<Cfg> nKrnl = deviceNKrnlPrototype_;
    nKrnl.numElements = (entry.get<real*>(inner_keys::Wp::Id::Dofs))->getSize();
    nKrnl.Qext = const_cast<const real**>(
        (entry.get<real*>(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr());
    nKrnl.Q = (entry.get<real*>(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
    nKrnl.Qane = (entry.get<real*>(inner_keys::Wp::Id::DofsAne))->getDeviceDataPtr();
    nKrnl.w = const_cast<const real**>(
        entry.get<real*>(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
    nKrnl.extraOffset_w = SEISSOL_OFFSET(LocalIntegrationData<Cfg>, specific.w);

    SEISSOL_OFFSET_ASSERT(LocalIntegrationData<Cfg>, specific.w);

    nKrnl.streamPtr = runtime.stream();

    nKrnl.execute();
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

#define SEISSOL_INSTANTIATE(Cfg) template class Neighbor<Cfg>;
SEISSOL_FOR_EACH_CONFIG_LINEARCKANELASTIC(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::kernels::solver::linearckanelastic
