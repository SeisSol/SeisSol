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
 * Local kernel of SeisSol.
 **/

#ifndef WAVEPROP_KERNEL_LOCAL_CKA_H_
#define WAVEPROP_KERNEL_LOCAL_CKA_H_

#include <memory>
#include "Common/configtensor.hpp"
#include "Physics/InitialField.h"
#include "Equations/datastructures.hpp"
#include <cassert>
#include <stdint.h>
#include <cstring>

#include <yateto.h>
#include "Equations/Local.hpp"

namespace seissol::waveprop::kernel::local {
template <typename Config>
class Local<Config,
            std::enable_if_t<Config::MaterialT::Solver ==
                             seissol::model::LocalSolver::CauchyKovalevskiAnelastic>> {
  using RealT = typename Config::RealT;

  protected:
  typename Yateto<Config>::Kernel::volumeExt m_volumeKernelPrototype;
  typename Yateto<Config>::Kernel::localFluxExt m_localFluxKernelPrototype;
  typename Yateto<Config>::Kernel::local m_localKernelPrototype;
  const std::vector<std::unique_ptr<physics::InitialField>>* initConds;

  public:
  virtual void setInitConds(decltype(initConds) initConds) { this->initConds = initConds; }

  void setHostGlobalData(GlobalData<Config> const* global) {
#ifndef NDEBUG
    for (unsigned stiffness = 0; stiffness < 3; ++stiffness) {
      assert(((uintptr_t)global->stiffnessMatrices(stiffness)) % Alignment == 0);
    }
    for (unsigned flux = 0; flux < 4; ++flux) {
      assert(((uintptr_t)global->localChangeOfBasisMatricesTransposed(flux)) % Alignment == 0);
      assert(((uintptr_t)global->changeOfBasisMatrices(flux)) % Alignment == 0);
    }
#endif

    m_volumeKernelPrototype.kDivM = global->stiffnessMatrices;
    m_localFluxKernelPrototype.rDivM = global->changeOfBasisMatrices;
    m_localFluxKernelPrototype.fMrT = global->localChangeOfBasisMatricesTransposed;
    m_localKernelPrototype.selectEla = Yateto<Config>::Init::selectEla::Values;
    m_localKernelPrototype.selectAne = Yateto<Config>::Init::selectAne::Values;
  }

  void setGlobalData(const CompoundGlobalData<Config>& global) { setHostGlobalData(global.onHost); }

  void computeIntegral(RealT i_timeIntegratedDegreesOfFreedom[Yateto<Config>::Tensor::I::size()],
                       LocalData<Config>& data,
                       LocalTmp<Config>& tmp,
                       // TODO(Lukas) Nullable cause miniseissol. Maybe fix?
                       const CellMaterialData* materialData,
                       CellBoundaryMapping const (*cellBoundaryMapping)[4],
                       double time,
                       double timeStepWidth) {
    // assert alignments
#ifndef NDEBUG
    assert(((uintptr_t)i_timeIntegratedDegreesOfFreedom) % Alignment == 0);
    assert(((uintptr_t)tmp.timeIntegratedAne) % Alignment == 0);
    assert(((uintptr_t)data.dofs) % Alignment == 0);
#endif

    RealT Qext[Yateto<Config>::Tensor::Qext::size()] __attribute__((aligned(Alignment)));

    typename Yateto<Config>::Kernel::volumeExt volKrnl = m_volumeKernelPrototype;
    volKrnl.Qext = Qext;
    volKrnl.I = i_timeIntegratedDegreesOfFreedom;
    for (unsigned i = 0; i < yateto::numFamilyMembers<typename Yateto<Config>::Tensor::star>();
         ++i) {
      volKrnl.star(i) = data.localIntegration.starMatrices[i];
    }

    typename Yateto<Config>::Kernel::localFluxExt lfKrnl = m_localFluxKernelPrototype;
    lfKrnl.Qext = Qext;
    lfKrnl.I = i_timeIntegratedDegreesOfFreedom;
    lfKrnl._prefetch.I = i_timeIntegratedDegreesOfFreedom + Yateto<Config>::Tensor::I::size();
    lfKrnl._prefetch.Q = data.dofs + Yateto<Config>::Tensor::Q::size();

    volKrnl.execute();

    for (unsigned int face = 0; face < 4; ++face) {
      // no element local contribution in the case of dynamic rupture boundary conditions
      if (data.cellInformation.faceTypes[face] != FaceType::dynamicRupture) {
        lfKrnl.AplusT = data.localIntegration.nApNm1[face];
        lfKrnl.execute(face);
      }
    }

    typename Yateto<Config>::Kernel::local lKrnl = m_localKernelPrototype;
    lKrnl.E = data.localIntegration.specific.E;
    lKrnl.Iane = tmp.timeIntegratedAne;
    lKrnl.Q = data.dofs;
    lKrnl.Qane = data.dofsAne;
    lKrnl.Qext = Qext;
    lKrnl.W = data.localIntegration.specific.W;
    lKrnl.w = data.localIntegration.specific.w;

    lKrnl.execute();
  }

  void flopsIntegral(FaceType const i_faceTypes[4],
                     unsigned int& o_nonZeroFlops,
                     unsigned int& o_hardwareFlops) {
    o_nonZeroFlops = seissol::Yateto<Config>::Kernel::volumeExt::NonZeroFlops;
    o_hardwareFlops = seissol::Yateto<Config>::Kernel::volumeExt::HardwareFlops;

    for (unsigned int face = 0; face < 4; ++face) {
      if (i_faceTypes[face] != FaceType::dynamicRupture) {
        o_nonZeroFlops += seissol::Yateto<Config>::Kernel::localFluxExt::nonZeroFlops(face);
        o_hardwareFlops += seissol::Yateto<Config>::Kernel::localFluxExt::hardwareFlops(face);
      }
    }

    o_nonZeroFlops += seissol::Yateto<Config>::Kernel::local::NonZeroFlops;
    o_hardwareFlops += seissol::Yateto<Config>::Kernel::local::HardwareFlops;
  }

  unsigned bytesIntegral() {
    unsigned reals = 0;

    // star matrices load
    reals += yateto::computeFamilySize<typename Yateto<Config>::Tensor::star>() +
             Yateto<Config>::Tensor::w::size() + Yateto<Config>::Tensor::W::size() +
             Yateto<Config>::Tensor::E::size();
    // flux solvers
    reals += 4 * Yateto<Config>::Tensor::AplusT::size();

    // DOFs write
    reals += Yateto<Config>::Tensor::Q::size() + Yateto<Config>::Tensor::Qane::size();

    return reals * sizeof(RealT);
  }
};
} // namespace seissol::waveprop::kernel::local

#endif
