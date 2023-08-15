/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de,
 *http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013-2014, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Neighbor kernel of SeisSol.
 **/

#ifndef WAVEPROP_KERNEL_NEIGHBOR_CKA_H_
#define WAVEPROP_KERNEL_NEIGHBOR_CKA_H_

#include <cassert>
#include <stdint.h>
#include <cstddef>
#include <cstring>

#include "Common/configtensor.hpp"
#include "Equations/datastructures.hpp"
#include "Equations/Neighbor.hpp"

namespace seissol::waveprop::kernel::neighbor {
template <typename Config>
class Neighbor<Config,
               std::enable_if_t<Config::MaterialT::Solver ==
                                seissol::model::LocalSolver::CauchyKovalevskiAnelastic>> {
  using RealT = typename Config::RealT;

  protected:
  typename Yateto<Config>::Kernel::neighbourFluxExt m_nfKrnlPrototype;
  typename Yateto<Config>::Kernel::neighbour m_nKrnlPrototype;
  typename Yateto<Config>::Kernel::dynamicRupture::nodalFlux m_drKrnlPrototype;

  public:
  void setHostGlobalData(GlobalData<Config> const* global) {
#ifndef NDEBUG
    for (int l_neighbor = 0; l_neighbor < 4; ++l_neighbor) {
      assert(((uintptr_t)global->changeOfBasisMatrices(l_neighbor)) % Alignment == 0);
      assert(((uintptr_t)global->localChangeOfBasisMatricesTransposed(l_neighbor)) % Alignment ==
             0);
      assert(((uintptr_t)global->neighbourChangeOfBasisMatricesTransposed(l_neighbor)) %
                 Alignment ==
             0);
    }

    for (int h = 0; h < 3; ++h) {
      assert(((uintptr_t)global->neighbourFluxMatrices(h)) % Alignment == 0);
    }

    for (int i = 0; i < 4; ++i) {
      for (int h = 0; h < 3; ++h) {
        assert(((uintptr_t)global->nodalFluxMatrices(i, h)) % Alignment == 0);
      }
    }
#endif
    m_nfKrnlPrototype.rDivM = global->changeOfBasisMatrices;
    m_nfKrnlPrototype.rT = global->neighbourChangeOfBasisMatricesTransposed;
    m_nfKrnlPrototype.fP = global->neighbourFluxMatrices;
    m_drKrnlPrototype.V3mTo2nTWDivM = global->nodalFluxMatrices;
    m_nKrnlPrototype.selectEla = Yateto<Config>::Init::selectEla::Values;
    m_nKrnlPrototype.selectAne = Yateto<Config>::Init::selectAne::Values;
  }

  void setGlobalData(const CompoundGlobalData<Config>& global) { setHostGlobalData(global.onHost); }

  void computeNeighborsIntegral(NeighborData<Config>& data,
                                CellDRMapping const (&cellDrMapping)[4],
                                RealT* i_timeIntegrated[4],
                                RealT* faceNeighbors_prefetch[4]) {
#ifndef NDEBUG
    for (int l_neighbor = 0; l_neighbor < 4; ++l_neighbor) {
      // alignment of the time integrated dofs
      if (data.cellInformation.faceTypes[l_neighbor] != FaceType::outflow &&
          data.cellInformation.faceTypes[l_neighbor] !=
              FaceType::dynamicRupture) { // no alignment for outflow and DR boundaries required
        assert(((uintptr_t)i_timeIntegrated[l_neighbor]) % Alignment == 0);
      }
    }
#endif

    // alignment of the degrees of freedom
    assert(((uintptr_t)data.dofs) % Alignment == 0);

    RealT Qext[Yateto<Config>::Tensor::Qext::size()] __attribute__((aligned(PAGESIZE_STACK))) = {};

    typename Yateto<Config>::Kernel::neighbourFluxExt nfKrnl = m_nfKrnlPrototype;
    nfKrnl.Qext = Qext;

    // iterate over faces
    for (unsigned int l_face = 0; l_face < 4; l_face++) {
      // no neighboring cell contribution in the case of absorbing and dynamic rupture boundary
      // conditions
      if (data.cellInformation.faceTypes[l_face] != FaceType::outflow &&
          data.cellInformation.faceTypes[l_face] != FaceType::dynamicRupture) {
        // compute the neighboring elements flux matrix id.
        if (data.cellInformation.faceTypes[l_face] != FaceType::freeSurface) {
          assert(data.cellInformation.faceRelations[l_face][0] < 4 &&
                 data.cellInformation.faceRelations[l_face][1] < 3);

          nfKrnl.I = i_timeIntegrated[l_face];
          nfKrnl.AminusT = data.neighboringIntegration.nAmNm1[l_face];
          nfKrnl._prefetch.I = faceNeighbors_prefetch[l_face];
          nfKrnl.execute(data.cellInformation.faceRelations[l_face][1],
                         data.cellInformation.faceRelations[l_face][0],
                         l_face);
        }
      } else if (data.cellInformation.faceTypes[l_face] == FaceType::dynamicRupture) {
        assert(((uintptr_t)cellDrMapping[l_face].godunov) % Alignment == 0);

        typename Yateto<Config>::Kernel::dynamicRupture::nodalFlux drKrnl = m_drKrnlPrototype;
        drKrnl.fluxSolver = cellDrMapping[l_face].fluxSolver;
        drKrnl.QInterpolated = cellDrMapping[l_face].godunov;
        drKrnl.Qext = Qext;
        drKrnl._prefetch.I = faceNeighbors_prefetch[l_face];
        drKrnl.execute(cellDrMapping[l_face].side, cellDrMapping[l_face].faceRelation);
      }
    }

    typename Yateto<Config>::Kernel::neighbour nKrnl = m_nKrnlPrototype;
    nKrnl.Qext = Qext;
    nKrnl.Q = data.dofs;
    nKrnl.Qane = data.dofsAne;
    nKrnl.w = data.neighboringIntegration.specific.w;

    nKrnl.execute();
  }

  void flopsNeighborsIntegral(const FaceType i_faceTypes[4],
                              const int i_neighboringIndices[4][2],
                              CellDRMapping const (&cellDrMapping)[4],
                              unsigned int& o_nonZeroFlops,
                              unsigned int& o_hardwareFlops,
                              long long& o_drNonZeroFlops,
                              long long& o_drHardwareFlops) {
    // reset flops
    o_nonZeroFlops = 0;
    o_hardwareFlops = 0;
    o_drNonZeroFlops = 0;
    o_drHardwareFlops = 0;

    for (unsigned int l_face = 0; l_face < 4; l_face++) {
      // no neighboring cell contribution in the case of absorbing and dynamic rupture boundary
      // conditions
      if (i_faceTypes[l_face] != FaceType::outflow &&
          i_faceTypes[l_face] != FaceType::dynamicRupture) {
        // compute the neighboring elements flux matrix id.
        if (i_faceTypes[l_face] != FaceType::freeSurface) {
          assert(i_neighboringIndices[l_face][0] < 4 && i_neighboringIndices[l_face][1] < 3);

          o_nonZeroFlops += seissol::Yateto<Config>::Kernel::neighbourFluxExt::nonZeroFlops(
              i_neighboringIndices[l_face][1], i_neighboringIndices[l_face][0], l_face);
          o_hardwareFlops += seissol::Yateto<Config>::Kernel::neighbourFluxExt::hardwareFlops(
              i_neighboringIndices[l_face][1], i_neighboringIndices[l_face][0], l_face);
        } else { // fall back to local matrices in case of free surface boundary conditions
          o_nonZeroFlops += seissol::Yateto<Config>::Kernel::localFluxExt::nonZeroFlops(l_face);
          o_hardwareFlops += seissol::Yateto<Config>::Kernel::localFluxExt::hardwareFlops(l_face);
        }
      } else if (i_faceTypes[l_face] == FaceType::dynamicRupture) {
        o_drNonZeroFlops += Yateto<Config>::Kernel::dynamicRupture::nodalFlux::nonZeroFlops(
            cellDrMapping[l_face].side, cellDrMapping[l_face].faceRelation);
        o_drHardwareFlops += Yateto<Config>::Kernel::dynamicRupture::nodalFlux::hardwareFlops(
            cellDrMapping[l_face].side, cellDrMapping[l_face].faceRelation);
      }
    }

    o_nonZeroFlops += Yateto<Config>::Kernel::neighbour::NonZeroFlops;
    o_hardwareFlops += Yateto<Config>::Kernel::neighbour::HardwareFlops;
  }

  unsigned bytesNeighborsIntegral() {
    unsigned reals = 0;

    // 4 * tElasticDOFS load, DOFs load, DOFs write
    reals += 4 * Yateto<Config>::Tensor::I::size() + 2 * Yateto<Config>::Tensor::Q::size() +
             2 * Yateto<Config>::Tensor::Qane::size();
    // flux solvers load
    reals += 4 * Yateto<Config>::Tensor::AminusT::size() + Yateto<Config>::Tensor::w::size();

    return reals * sizeof(RealT);
  }
};
} // namespace seissol::waveprop::kernel::neighbor

#endif
