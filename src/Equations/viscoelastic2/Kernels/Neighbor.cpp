// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#include "Kernels/Neighbor.h"

#include <cassert>
#include <cstddef>
#include <cstring>
#include <stdint.h>

#include "generated_code/init.h"

namespace seissol::kernels {

void Neighbor::setHostGlobalData(const GlobalData* global) {
#ifndef NDEBUG
  for (int neighbor = 0; neighbor < 4; ++neighbor) {
    assert((reinterpret_cast<uintptr_t>(global->changeOfBasisMatrices(neighbor))) % Alignment == 0);
    assert((reinterpret_cast<uintptr_t>(global->localChangeOfBasisMatricesTransposed(neighbor))) %
               Alignment ==
           0);
    assert(
        (reinterpret_cast<uintptr_t>(global->neighbourChangeOfBasisMatricesTransposed(neighbor))) %
            Alignment ==
        0);
  }

  for (int h = 0; h < 3; ++h) {
    assert((reinterpret_cast<uintptr_t>(global->neighbourFluxMatrices(h))) % Alignment == 0);
  }

  for (int i = 0; i < 4; ++i) {
    for (int h = 0; h < 3; ++h) {
      assert((reinterpret_cast<uintptr_t>(global->nodalFluxMatrices(i, h))) % Alignment == 0);
    }
  }
#endif
  m_nfKrnlPrototype.rDivM = global->changeOfBasisMatrices;
  m_nfKrnlPrototype.rT = global->neighbourChangeOfBasisMatricesTransposed;
  m_nfKrnlPrototype.fP = global->neighbourFluxMatrices;
  m_drKrnlPrototype.V3mTo2nTWDivM = global->nodalFluxMatrices;
  m_nKrnlPrototype.selectEla = init::selectEla::Values;
  m_nKrnlPrototype.selectAne = init::selectAne::Values;
}

void Neighbor::setGlobalData(const CompoundGlobalData& global) { setHostGlobalData(global.onHost); }

void Neighbor::computeNeighborsIntegral(NeighborData& data,
                                        const CellDRMapping (&cellDrMapping)[4],
                                        real* timeIntegrated[4],
                                        real* faceNeighbors_prefetch[4]) {
#ifndef NDEBUG
  for (int neighbor = 0; neighbor < 4; ++neighbor) {
    // alignment of the time integrated dofs
    if (data.cellInformation().faceTypes[neighbor] != FaceType::Outflow &&
        data.cellInformation().faceTypes[neighbor] !=
            FaceType::DynamicRupture) { // no alignment for outflow and DR boundaries required
      assert((reinterpret_cast<uintptr_t>(timeIntegrated[neighbor])) % Alignment == 0);
    }
  }
#endif

  // alignment of the degrees of freedom
  assert((reinterpret_cast<uintptr_t>(data.dofs())) % Alignment == 0);

  alignas(PagesizeStack) real Qext[tensor::Qext::size()] = {};

  kernel::neighbourFluxExt nfKrnl = m_nfKrnlPrototype;
  nfKrnl.Qext = Qext;

  // iterate over faces
  for (unsigned int face = 0; face < 4; face++) {
    // no neighboring cell contribution in the case of absorbing and dynamic rupture boundary
    // conditions
    if (data.cellInformation().faceTypes[face] != FaceType::Outflow &&
        data.cellInformation().faceTypes[face] != FaceType::DynamicRupture) {
      // compute the neighboring elements flux matrix id.
      if (data.cellInformation().faceTypes[face] != FaceType::FreeSurface) {
        assert(data.cellInformation().faceRelations[face][0] < 4 &&
               data.cellInformation().faceRelations[face][1] < 3);

        nfKrnl.I = timeIntegrated[face];
        nfKrnl.AminusT = data.neighboringIntegration().nAmNm1[face];
        nfKrnl._prefetch.I = faceNeighbors_prefetch[face];
        nfKrnl.execute(data.cellInformation().faceRelations[face][1],
                       data.cellInformation().faceRelations[face][0],
                       face);
      }
    } else if (data.cellInformation().faceTypes[face] == FaceType::DynamicRupture) {
      assert((reinterpret_cast<uintptr_t>(cellDrMapping[face].godunov)) % Alignment == 0);

      dynamicRupture::kernel::nodalFlux drKrnl = m_drKrnlPrototype;
      drKrnl.fluxSolver = cellDrMapping[face].fluxSolver;
      drKrnl.QInterpolated = cellDrMapping[face].godunov;
      drKrnl.Qext = Qext;
      drKrnl._prefetch.I = faceNeighbors_prefetch[face];
      drKrnl.execute(cellDrMapping[face].side, cellDrMapping[face].faceRelation);
    }
  }

  kernel::neighbour nKrnl = m_nKrnlPrototype;
  nKrnl.Qext = Qext;
  nKrnl.Q = data.dofs();
  nKrnl.Qane = data.dofsAne();
  nKrnl.w = data.neighboringIntegration().specific.w;

  nKrnl.execute();
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
    // no neighboring cell contribution in the case of absorbing and dynamic rupture boundary
    // conditions
    if (faceTypes[face] != FaceType::Outflow && faceTypes[face] != FaceType::DynamicRupture) {
      // compute the neighboring elements flux matrix id.
      if (faceTypes[face] != FaceType::FreeSurface) {
        assert(neighboringIndices[face][0] < 4 && neighboringIndices[face][1] < 3);

        nonZeroFlops += seissol::kernel::neighbourFluxExt::nonZeroFlops(
            neighboringIndices[face][1], neighboringIndices[face][0], face);
        hardwareFlops += seissol::kernel::neighbourFluxExt::hardwareFlops(
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

  nonZeroFlops += kernel::neighbour::NonZeroFlops;
  hardwareFlops += kernel::neighbour::HardwareFlops;
}

unsigned Neighbor::bytesNeighborsIntegral() {
  unsigned reals = 0;

  // 4 * tElasticDOFS load, DOFs load, DOFs write
  reals += 4 * tensor::I::size() + 2 * tensor::Q::size() + 2 * tensor::Qane::size();
  // flux solvers load
  reals += 4 * tensor::AminusT::size() + tensor::w::size();

  return reals * sizeof(real);
}

} // namespace seissol::kernels
