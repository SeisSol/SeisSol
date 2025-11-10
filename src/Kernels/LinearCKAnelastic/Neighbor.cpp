// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#include "NeighborBase.h"

#include <Alignment.h>
#include <Common/Constants.h>
#include <DataTypes/ConditionalTable.h>
#include <GeneratedCode/metagen/kernel.h>
#include <GeneratedCode/metagen/tensor.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/Typedefs.h>
#include <Memory/Descriptor/LTS.h>
#include <Parallel/Runtime/Stream.h>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdint.h>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#endif

#include "GeneratedCode/init.h"

#include "Memory/GlobalData.h"

namespace seissol::kernels::solver::linearckanelastic {

template <typename Cfg>
void Neighbor<Cfg>::setGlobalData(const GlobalData& global) {
#ifndef NDEBUG
  for (std::size_t neighbor = 0; neighbor < Cell::NumFaces; ++neighbor) {
    assert((reinterpret_cast<uintptr_t>(global.get<Cfg>().changeOfBasisMatrices(neighbor))) %
               Alignment ==
           0);
    assert((reinterpret_cast<uintptr_t>(
               global.get<Cfg>().localChangeOfBasisMatricesTransposed(neighbor))) %
               Alignment ==
           0);
    assert((reinterpret_cast<uintptr_t>(
               global.get<Cfg>().neighborChangeOfBasisMatricesTransposed(neighbor))) %
               Alignment ==
           0);
  }

  for (int h = 0; h < 3; ++h) {
    assert((reinterpret_cast<uintptr_t>(global.get<Cfg>().neighborFluxMatrices(h))) % Alignment ==
           0);
  }

  for (std::size_t i = 0; i < Cell::NumFaces; ++i) {
    for (int h = 0; h < 3; ++h) {
      assert((reinterpret_cast<uintptr_t>(global.get<Cfg>().nodalFluxMatrices(i, h))) % Alignment ==
             0);
    }
  }
#endif
  m_nfKrnlPrototype.rDivM = global.get<Cfg>().changeOfBasisMatrices;
  m_nfKrnlPrototype.rT = global.get<Cfg>().neighborChangeOfBasisMatricesTransposed;
  m_nfKrnlPrototype.fP = global.get<Cfg>().neighborFluxMatrices;
  m_drKrnlPrototype.V3mTo2nTWDivM = global.get<Cfg>().nodalFluxMatrices;
  m_nKrnlPrototype.selectEla = init::selectEla<Cfg>::Values;
  m_nKrnlPrototype.selectAne = init::selectAne<Cfg>::Values;

#ifdef ACL_DEVICE
#ifdef USE_PREMULTIPLY_FLUX
  deviceNfKrnlPrototype.minusFluxMatrices = global.get<Cfg, Executor::Device>().minusFluxMatrices;
#else
  deviceNfKrnlPrototype.rDivM = global.get<Cfg, Executor::Device>().changeOfBasisMatrices;
  deviceNfKrnlPrototype.rT =
      global.get<Cfg, Executor::Device>().neighborChangeOfBasisMatricesTransposed;
  deviceNfKrnlPrototype.fP = global.get<Cfg, Executor::Device>().neighborFluxMatrices;
#endif
  deviceDrKrnlPrototype.V3mTo2nTWDivM = global.get<Cfg, Executor::Device>().nodalFluxMatrices;
  deviceNKrnlPrototype.selectEla = global.get<Cfg, Executor::Device>().selectEla;
  deviceNKrnlPrototype.selectAne = global.get<Cfg, Executor::Device>().selectAne;
#endif
}

template <typename Cfg>
void Neighbor<Cfg>::computeNeighborsIntegral(LTS::Ref<Cfg>& data,
                                             const CellDRMapping<Cfg> (&cellDrMapping)[4],
                                             real* timeIntegrated[4],
                                             real* faceNeighborsPrefetch[4]) {
#ifndef NDEBUG
  for (std::size_t neighbor = 0; neighbor < Cell::NumFaces; ++neighbor) {
    // alignment of the time integrated dofs
    if (data.template get<LTS::CellInformation>().faceTypes[neighbor] != FaceType::Outflow &&
        data.template get<LTS::CellInformation>().faceTypes[neighbor] !=
            FaceType::DynamicRupture) { // no alignment for outflow and DR boundaries required
      assert((reinterpret_cast<uintptr_t>(timeIntegrated[neighbor])) % Alignment == 0);
    }
  }
#endif

  // alignment of the degrees of freedom
  assert((reinterpret_cast<uintptr_t>(data.template get<LTS::Dofs>())) % Alignment == 0);

  alignas(PagesizeStack) real qext[tensor::Qext<Cfg>::size()] = {};

  kernel::neighborFluxExt<Cfg> nfKrnl = m_nfKrnlPrototype;
  nfKrnl.Qext = qext;

  // iterate over faces
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // no neighboring cell contribution in the case of absorbing and dynamic rupture boundary
    // conditions
    if (data.template get<LTS::CellInformation>().faceTypes[face] != FaceType::Outflow &&
        data.template get<LTS::CellInformation>().faceTypes[face] != FaceType::DynamicRupture) {
      // compute the neighboring elements flux matrix id.
      if (data.template get<LTS::CellInformation>().faceTypes[face] != FaceType::FreeSurface) {
        assert(data.template get<LTS::CellInformation>().faceRelations[face][0] < Cell::NumFaces &&
               data.template get<LTS::CellInformation>().faceRelations[face][1] < 3);

        nfKrnl.I = timeIntegrated[face];
        nfKrnl.AminusT = data.template get<LTS::NeighboringIntegration>().nAmNm1[face];
        nfKrnl._prefetch.I = faceNeighborsPrefetch[face];
        nfKrnl.execute(data.template get<LTS::CellInformation>().faceRelations[face][1],
                       data.template get<LTS::CellInformation>().faceRelations[face][0],
                       face);
      }
    } else if (data.template get<LTS::CellInformation>().faceTypes[face] ==
               FaceType::DynamicRupture) {
      assert((reinterpret_cast<uintptr_t>(cellDrMapping[face].godunov)) % Alignment == 0);

      dynamicRupture::kernel::nodalFlux<Cfg> drKrnl = m_drKrnlPrototype;
      drKrnl.fluxSolver = cellDrMapping[face].fluxSolver;
      drKrnl.QInterpolated = cellDrMapping[face].godunov;
      drKrnl.Qext = qext;
      drKrnl._prefetch.I = faceNeighborsPrefetch[face];
      drKrnl.execute(cellDrMapping[face].side, cellDrMapping[face].faceRelation);
    }
  }

  kernel::neighbor<Cfg> nKrnl = m_nKrnlPrototype;
  nKrnl.Qext = qext;
  nKrnl.Q = data.template get<LTS::Dofs>();
  nKrnl.Qane = data.template get<LTS::DofsAne>();
  nKrnl.w = data.template get<LTS::NeighboringIntegration>().specific.w;

  nKrnl.execute();
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
    // no neighboring cell contribution in the case of absorbing and dynamic rupture boundary
    // conditions
    if (faceTypes[face] != FaceType::Outflow && faceTypes[face] != FaceType::DynamicRupture) {
      // compute the neighboring elements flux matrix id.
      if (faceTypes[face] != FaceType::FreeSurface) {
        assert(neighboringIndices[face][0] < Cell::NumFaces && neighboringIndices[face][1] < 3);

        nonZeroFlops += seissol::kernel::neighborFluxExt<Cfg>::nonZeroFlops(
            neighboringIndices[face][1], neighboringIndices[face][0], face);
        hardwareFlops += seissol::kernel::neighborFluxExt<Cfg>::hardwareFlops(
            neighboringIndices[face][1], neighboringIndices[face][0], face);
      } else { // fall back to local matrices in case of free surface boundary conditions
        nonZeroFlops += seissol::kernel::localFluxExt<Cfg>::nonZeroFlops(face);
        hardwareFlops += seissol::kernel::localFluxExt<Cfg>::hardwareFlops(face);
      }
    } else if (faceTypes[face] == FaceType::DynamicRupture) {
      drNonZeroFlops += dynamicRupture::kernel::nodalFlux<Cfg>::nonZeroFlops(
          cellDrMapping[face].side, cellDrMapping[face].faceRelation);
      drHardwareFlops += dynamicRupture::kernel::nodalFlux<Cfg>::hardwareFlops(
          cellDrMapping[face].side, cellDrMapping[face].faceRelation);
    }
  }

  nonZeroFlops += kernel::neighbor<Cfg>::NonZeroFlops;
  hardwareFlops += kernel::neighbor<Cfg>::HardwareFlops;
}

template <typename Cfg>
std::uint64_t Neighbor<Cfg>::bytesNeighborsIntegral() {
  std::uint64_t reals = 0;

  // 4 * tElasticDOFS load, DOFs load, DOFs write
  reals += 4 * tensor::I<Cfg>::size() + 2 * tensor::Q<Cfg>::size() + 2 * tensor::Qane<Cfg>::size();
  // flux solvers load
  reals += 4 * tensor::AminusT<Cfg>::size() + tensor::w<Cfg>::size();

  return reals * sizeof(real);
}

template <typename Cfg>
void Neighbor<Cfg>::computeBatchedNeighborsIntegral(
    ConditionalPointersToRealsTable& table, seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  kernel::gpu_neighborFluxExt<Cfg> neighFluxKrnl = deviceNfKrnlPrototype;
  dynamicRupture::kernel::gpu_nodalFlux<Cfg> drKrnl = deviceDrKrnlPrototype;

  {
    ConditionalKey key(KernelNames::Time || KernelNames::Volume);
    if (table.find(key) != table.end()) {
      auto& entry = table[key];
      device.algorithms.setToValue((entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr(),
                                   static_cast<real>(0.0),
                                   tensor::Qext<Cfg>::Size,
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
            unsigned faceRelation = i;

            ConditionalKey key(*KernelNames::NeighborFlux,
                               (FaceKinds::Regular || FaceKinds::Periodic),
                               face,
                               faceRelation);

            if (table.find(key) != table.end()) {
              auto& entry = table[key];

              const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
              neighFluxKrnl.numElements = numElements;

              neighFluxKrnl.Qext = (entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr();
              neighFluxKrnl.I = const_cast<const real**>(
                  (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());
              neighFluxKrnl.AminusT = const_cast<const real**>(
                  entry.get(inner_keys::Wp::Id::NeighborIntegrationData)->getDeviceDataPtr());
              neighFluxKrnl.extraOffset_AminusT =
                  SEISSOL_ARRAY_OFFSET(NeighboringIntegrationData, nAmNm1, face);

              neighFluxKrnl.streamPtr = stream;
              (neighFluxKrnl.*
               seissol::kernel::gpu_neighborFluxExt<Cfg>::ExecutePtrs[faceRelation])();
            }
          } else {
            // Dynamic Rupture
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
              drKrnl.Qext = (entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr();

              drKrnl.streamPtr = stream;
              (drKrnl.*
               seissol::dynamicRupture::kernel::gpu_nodalFlux<Cfg>::ExecutePtrs[faceRelation])();
            }
          }
        });
  }

  ConditionalKey key(KernelNames::Time || KernelNames::Volume);
  if (table.find(key) != table.end()) {
    auto& entry = table[key];
    kernel::gpu_neighbor<Cfg> nKrnl = deviceNKrnlPrototype;
    nKrnl.numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
    nKrnl.Qext =
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr());
    nKrnl.Q = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
    nKrnl.Qane = (entry.get(inner_keys::Wp::Id::DofsAne))->getDeviceDataPtr();
    nKrnl.w = const_cast<const real**>(
        entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
    nKrnl.extraOffset_w = SEISSOL_OFFSET(LocalIntegrationData, specific.w);

    nKrnl.streamPtr = runtime.stream();

    nKrnl.execute();
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

#define SEISSOL_CONFIGITER(cfg) template class Neighbor<cfg>;
#include "ConfigIncludeLinearCKAne.h"

} // namespace seissol::kernels::solver::linearckanelastic
