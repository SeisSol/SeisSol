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

#include <cassert>
#include <cstddef>
#include <cstring>
#include <stdint.h>

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#endif

namespace seissol::kernels::solver::linearckanelastic {

void Neighbor::setGlobalData(const CompoundGlobalData& global) {
#ifndef NDEBUG
  for (std::size_t neighbor = 0; neighbor < Cell::NumFaces; ++neighbor) {
    assert((reinterpret_cast<uintptr_t>(global.onHost->changeOfBasisMatrices(neighbor))) %
               Alignment ==
           0);
    assert((reinterpret_cast<uintptr_t>(
               global.onHost->localChangeOfBasisMatricesTransposed(neighbor))) %
               Alignment ==
           0);
    assert((reinterpret_cast<uintptr_t>(
               global.onHost->neighborChangeOfBasisMatricesTransposed(neighbor))) %
               Alignment ==
           0);
  }

  for (int h = 0; h < 3; ++h) {
    assert((reinterpret_cast<uintptr_t>(global.onHost->neighborFluxMatrices(h))) % Alignment == 0);
  }

  for (std::size_t i = 0; i < Cell::NumFaces; ++i) {
    for (int h = 0; h < 3; ++h) {
      assert((reinterpret_cast<uintptr_t>(global.onHost->nodalFluxMatrices(i, h))) % Alignment ==
             0);
    }
  }
#endif
  m_nfKrnlPrototype.rDivM = global.onHost->changeOfBasisMatrices;
  m_nfKrnlPrototype.rT = global.onHost->neighborChangeOfBasisMatricesTransposed;
  m_nfKrnlPrototype.fP = global.onHost->neighborFluxMatrices;
  m_drKrnlPrototype.V3mTo2nTWDivM = global.onHost->nodalFluxMatrices;

#ifdef ACL_DEVICE
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

void Neighbor::computeNeighborsIntegral(LTS::Ref& data,
                                        const CellDRMapping (&cellDrMapping)[4],
                                        real* timeIntegrated[4],
                                        real* faceNeighbors_prefetch[4]) {
#ifndef NDEBUG
  for (std::size_t neighbor = 0; neighbor < Cell::NumFaces; ++neighbor) {
    // alignment of the time integrated dofs
    if (data.get<LTS::CellInformation>().faceTypes[neighbor] != FaceType::Outflow &&
        data.get<LTS::CellInformation>().faceTypes[neighbor] !=
            FaceType::DynamicRupture) { // no alignment for outflow and DR boundaries required
      assert((reinterpret_cast<uintptr_t>(timeIntegrated[neighbor])) % Alignment == 0);
    }
  }
#endif

  // alignment of the degrees of freedom
  assert((reinterpret_cast<uintptr_t>(data.get<LTS::Dofs>())) % Alignment == 0);

  alignas(PagesizeStack) real Qext[tensor::Qext::size()] = {};

  kernel::neighborFluxExt nfKrnl = m_nfKrnlPrototype;
  nfKrnl.Qext = Qext;

  // iterate over faces
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // no neighboring cell contribution in the case of absorbing and dynamic rupture boundary
    // conditions
    if (data.get<LTS::CellInformation>().faceTypes[face] != FaceType::Outflow &&
        data.get<LTS::CellInformation>().faceTypes[face] != FaceType::DynamicRupture) {
      // compute the neighboring elements flux matrix id.
      if (data.get<LTS::CellInformation>().faceTypes[face] != FaceType::FreeSurface) {
        assert(data.get<LTS::CellInformation>().faceRelations[face][0] < Cell::NumFaces &&
               data.get<LTS::CellInformation>().faceRelations[face][1] < 3);

        nfKrnl.I = timeIntegrated[face];
        nfKrnl.AminusT = data.get<LTS::NeighboringIntegration>().nAmNm1[face];
        nfKrnl._prefetch.I = faceNeighbors_prefetch[face];
        nfKrnl.execute(data.get<LTS::CellInformation>().faceRelations[face][1],
                       data.get<LTS::CellInformation>().faceRelations[face][0],
                       face);
      }
    } else if (data.get<LTS::CellInformation>().faceTypes[face] == FaceType::DynamicRupture) {
      assert((reinterpret_cast<uintptr_t>(cellDrMapping[face].godunov)) % Alignment == 0);

      dynamicRupture::kernel::nodalFlux drKrnl = m_drKrnlPrototype;
      drKrnl.fluxSolver = cellDrMapping[face].fluxSolver;
      drKrnl.QInterpolated = cellDrMapping[face].godunov;
      drKrnl.Qext = Qext;
      drKrnl._prefetch.I = faceNeighbors_prefetch[face];
      drKrnl.execute(cellDrMapping[face].side, cellDrMapping[face].faceRelation);
    }
  }

  kernel::neighbor nKrnl = m_nKrnlPrototype;
  nKrnl.Qext = Qext;
  nKrnl.Q = data.get<LTS::Dofs>();
  nKrnl.Qane = data.get<LTS::DofsAne>();
  nKrnl.w = data.get<LTS::NeighboringIntegration>().specific.w;

  nKrnl.execute();
}

void Neighbor::flopsNeighborsIntegral(
    const std::array<FaceType, Cell::NumFaces>& faceTypes,
    const std::array<std::array<uint8_t, 2>, Cell::NumFaces>& neighboringIndices,
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
    // no neighboring cell contribution in the case of absorbing and dynamic rupture boundary
    // conditions
    if (faceTypes[face] != FaceType::Outflow && faceTypes[face] != FaceType::DynamicRupture) {
      // compute the neighboring elements flux matrix id.
      if (faceTypes[face] != FaceType::FreeSurface) {
        assert(neighboringIndices[face][0] < Cell::NumFaces && neighboringIndices[face][1] < 3);

        nonZeroFlops += seissol::kernel::neighborFluxExt::nonZeroFlops(
            neighboringIndices[face][1], neighboringIndices[face][0], face);
        hardwareFlops += seissol::kernel::neighborFluxExt::hardwareFlops(
            neighboringIndices[face][1], neighboringIndices[face][0], face);
      } else { // fall back to local matrices in case of free surface boundary conditions
        nonZeroFlops += seissol::kernel::localFluxExt::nonZeroFlops(face);
        hardwareFlops += seissol::kernel::localFluxExt::hardwareFlops(face);
      }
    } else if (faceTypes[face] == FaceType::DynamicRupture) {
      drNonZeroFlops += dynamicRupture::kernel::nodalFlux::nonZeroFlops(
          cellDrMapping[face].side, cellDrMapping[face].faceRelation);
      drHardwareFlops += dynamicRupture::kernel::nodalFlux::hardwareFlops(
          cellDrMapping[face].side, cellDrMapping[face].faceRelation);
    }
  }

  nonZeroFlops += kernel::neighbor::NonZeroFlops;
  hardwareFlops += kernel::neighbor::HardwareFlops;
}

std::uint64_t Neighbor::bytesNeighborsIntegral() {
  std::uint64_t reals = 0;

  // 4 * tElasticDOFS load, DOFs load, DOFs write
  reals += 4 * tensor::I::size() + 2 * tensor::Q::size() + 2 * tensor::Qane::size();
  // flux solvers load
  reals += 4 * tensor::AminusT::size() + tensor::w::size();

  return reals * sizeof(real);
}

void Neighbor::computeBatchedNeighborsIntegral(
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& table,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE

  using namespace seissol::recording;
  kernel::gpu_neighborFluxExt neighFluxKrnl = deviceNfKrnlPrototype;
  dynamicRupture::kernel::gpu_nodalFlux drKrnl = deviceDrKrnlPrototype;

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
              (neighFluxKrnl.*seissol::kernel::gpu_neighborFluxExt::ExecutePtrs[faceRelation])();
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
              (drKrnl.*seissol::dynamicRupture::kernel::gpu_nodalFlux::ExecutePtrs[faceRelation])();
            }
          }
        });
  }

  ConditionalKey key(KernelNames::Time || KernelNames::Volume);
  if (table.find(key) != table.end()) {
    auto& entry = table[key];
    kernel::gpu_neighbor nKrnl = deviceNKrnlPrototype;
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

} // namespace seissol::kernels::solver::linearckanelastic
