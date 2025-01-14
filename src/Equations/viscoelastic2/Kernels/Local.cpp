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

#include "Kernels/Local.h"

#include <cassert>
#include <cstring>
#include <stdint.h>

#include <yateto.h>

namespace seissol::kernels {

void Local::setHostGlobalData(const GlobalData* global) {
#ifndef NDEBUG
  for (unsigned stiffness = 0; stiffness < 3; ++stiffness) {
    assert((reinterpret_cast<uintptr_t>(global->stiffnessMatrices(stiffness))) % Alignment == 0);
  }
  for (unsigned flux = 0; flux < 4; ++flux) {
    assert((reinterpret_cast<uintptr_t>(global->localChangeOfBasisMatricesTransposed(flux))) %
               Alignment ==
           0);
    assert((reinterpret_cast<uintptr_t>(global->changeOfBasisMatrices(flux))) % Alignment == 0);
  }
#endif

  m_volumeKernelPrototype.kDivM = global->stiffnessMatrices;
  m_localFluxKernelPrototype.rDivM = global->changeOfBasisMatrices;
  m_localFluxKernelPrototype.fMrT = global->localChangeOfBasisMatricesTransposed;
  m_localKernelPrototype.selectEla = init::selectEla::Values;
  m_localKernelPrototype.selectAne = init::selectAne::Values;
}

void Local::setGlobalData(const CompoundGlobalData& global) { setHostGlobalData(global.onHost); }

void Local::computeIntegral(real timeIntegratedDegreesOfFreedom[tensor::I::size()],
                            LocalData& data,
                            LocalTmp& tmp,
                            // TODO(Lukas) Nullable cause miniseissol. Maybe fix?
                            const CellMaterialData* materialData,
                            const CellBoundaryMapping (*cellBoundaryMapping)[4],
                            double time,
                            double timeStepWidth) {
  // assert alignments
#ifndef NDEBUG
  assert((reinterpret_cast<uintptr_t>(timeIntegratedDegreesOfFreedom)) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(tmp.timeIntegratedAne)) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(data.dofs())) % Alignment == 0);
#endif

  alignas(Alignment) real Qext[tensor::Qext::size()];

  kernel::volumeExt volKrnl = m_volumeKernelPrototype;
  volKrnl.Qext = Qext;
  volKrnl.I = timeIntegratedDegreesOfFreedom;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    volKrnl.star(i) = data.localIntegration().starMatrices[i];
  }

  kernel::localFluxExt lfKrnl = m_localFluxKernelPrototype;
  lfKrnl.Qext = Qext;
  lfKrnl.I = timeIntegratedDegreesOfFreedom;
  lfKrnl._prefetch.I = timeIntegratedDegreesOfFreedom + tensor::I::size();
  lfKrnl._prefetch.Q = data.dofs() + tensor::Q::size();

  volKrnl.execute();

  for (unsigned int face = 0; face < 4; ++face) {
    // no element local contribution in the case of dynamic rupture boundary conditions
    if (data.cellInformation().faceTypes[face] != FaceType::DynamicRupture) {
      lfKrnl.AplusT = data.localIntegration().nApNm1[face];
      lfKrnl.execute(face);
    }
  }

  kernel::local lKrnl = m_localKernelPrototype;
  lKrnl.E = data.localIntegration().specific.E;
  lKrnl.Iane = tmp.timeIntegratedAne;
  lKrnl.Q = data.dofs();
  lKrnl.Qane = data.dofsAne();
  lKrnl.Qext = Qext;
  lKrnl.W = data.localIntegration().specific.W;
  lKrnl.w = data.localIntegration().specific.w;

  lKrnl.execute();
}

void Local::flopsIntegral(const FaceType faceTypes[4],
                          unsigned int& nonZeroFlops,
                          unsigned int& hardwareFlops) {
  nonZeroFlops = seissol::kernel::volumeExt::NonZeroFlops;
  hardwareFlops = seissol::kernel::volumeExt::HardwareFlops;

  for (unsigned int face = 0; face < 4; ++face) {
    if (faceTypes[face] != FaceType::DynamicRupture) {
      nonZeroFlops += seissol::kernel::localFluxExt::nonZeroFlops(face);
      hardwareFlops += seissol::kernel::localFluxExt::hardwareFlops(face);
    }
  }

  nonZeroFlops += seissol::kernel::local::NonZeroFlops;
  hardwareFlops += seissol::kernel::local::HardwareFlops;
}

unsigned Local::bytesIntegral() {
  unsigned reals = 0;

  // star matrices load
  reals += yateto::computeFamilySize<tensor::star>() + tensor::w::size() + tensor::W::size() +
           tensor::E::size();
  // flux solvers
  reals += 4 * tensor::AplusT::size();

  // DOFs write
  reals += tensor::Q::size() + tensor::Qane::size();

  return reals * sizeof(real);
}

} // namespace seissol::kernels
