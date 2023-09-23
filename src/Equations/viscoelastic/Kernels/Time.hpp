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
 * Copyright (c) 2013-2015, SeisSol Group
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
 * Time kernel of SeisSol.
 **/

#ifndef WAVEPROP_KERNEL_TIME_CKA_H_
#define WAVEPROP_KERNEL_TIME_CKA_H_

#include <cstring>
#include <cassert>
#include <stdint.h>

#include <yateto.h>

#include <Kernels/denseMatrixOps.hpp>
#include "Equations/datastructures.hpp"
#include "Equations/Time.hpp"
#include "Common/configtensor.hpp"

namespace seissol::waveprop::kernel::time {
template <typename Config>
class Time<Config,
           std::enable_if_t<Config::MaterialT::Solver ==
                            seissol::model::LocalSolver::CauchyKovalevskiAnelastic>> {
  using RealT = typename Config::RealT;

  protected:
  typename Yateto<Config>::Kernel::derivative m_krnlPrototype;

  public:
  void setHostGlobalData(GlobalData<Config> const* global) {
    assert(((uintptr_t)global->stiffnessMatricesTransposed(0)) % Alignment == 0);
    assert(((uintptr_t)global->stiffnessMatricesTransposed(1)) % Alignment == 0);
    assert(((uintptr_t)global->stiffnessMatricesTransposed(2)) % Alignment == 0);

    m_krnlPrototype.kDivMT = global->stiffnessMatricesTransposed;
    m_krnlPrototype.selectAne = Yateto<Config>::Init::selectAne::Values;
    m_krnlPrototype.selectEla = Yateto<Config>::Init::selectEla::Values;
  }

  void setGlobalData(const CompoundGlobalData<Config>& global) { setHostGlobalData(global.onHost); }

  void computeAder(double i_timeStepWidth,
                   LocalData<Config>& data,
                   LocalTmp<Config>& tmp,
                   RealT o_timeIntegrated[ConfigConstants<Config>::TensorSizeI],
                   RealT* o_timeDerivatives,
                   bool updateDisplacement) {
    /*
     * assert alignments.
     */
    assert(((uintptr_t)data.dofs) % Alignment == 0);
    assert(((uintptr_t)o_timeIntegrated) % Alignment == 0);
    assert(((uintptr_t)o_timeDerivatives) % Alignment == 0 || o_timeDerivatives == NULL);

    /*
     * compute ADER scheme.
     */
    // temporary result
    RealT temporaryBuffer[2][Yateto<Config>::Tensor::dQ::size(0)]
        __attribute__((aligned(PAGESIZE_STACK)));
    RealT temporaryBufferExt[2][Yateto<Config>::Tensor::dQext::size(1)]
        __attribute__((aligned(PAGESIZE_STACK)));
    RealT temporaryBufferAne[2][Yateto<Config>::Tensor::dQane::size(0)]
        __attribute__((aligned(PAGESIZE_STACK)));

    typename Yateto<Config>::Kernel::derivative krnl = m_krnlPrototype;
    typename Yateto<Config>::Kernel::derivativeTaylorExpansion intKrnl;

    krnl.dQ(0) = const_cast<RealT*>(data.dofs);
    intKrnl.dQ(0) = data.dofs;
    if (o_timeDerivatives != nullptr) {
      streamstore(Yateto<Config>::Tensor::dQ::size(0), data.dofs, o_timeDerivatives);
      RealT* derOut = o_timeDerivatives;
      for (unsigned i = 1; i < yateto::numFamilyMembers<typename Yateto<Config>::Tensor::dQ>();
           ++i) {
        derOut += Yateto<Config>::Tensor::dQ::size(i - 1);
        krnl.dQ(i) = derOut;
        intKrnl.dQ(i) = derOut;
      }
    } else {
      for (unsigned i = 1; i < yateto::numFamilyMembers<typename Yateto<Config>::Tensor::dQ>();
           ++i) {
        krnl.dQ(i) = temporaryBuffer[i % 2];
        intKrnl.dQ(i) = temporaryBuffer[i % 2];
      }
    }

    krnl.dQane(0) = const_cast<RealT*>(data.dofsAne);
    intKrnl.dQane(0) = const_cast<RealT*>(data.dofsAne);
    for (unsigned i = 1; i < yateto::numFamilyMembers<typename Yateto<Config>::Tensor::dQ>(); ++i) {
      krnl.dQane(i) = temporaryBufferAne[i % 2];
      intKrnl.dQane(i) = temporaryBufferAne[i % 2];
      krnl.dQext(i) = temporaryBufferExt[i % 2];
    }

    intKrnl.I = o_timeIntegrated;
    intKrnl.Iane = tmp.timeIntegratedAne;

    for (unsigned i = 0; i < yateto::numFamilyMembers<typename Yateto<Config>::Tensor::star>();
         ++i) {
      krnl.star(i) = data.localIntegration.starMatrices[i];
    }
    krnl.w = data.localIntegration.specific.w;
    krnl.W = data.localIntegration.specific.W;
    krnl.E = data.localIntegration.specific.E;

    // powers in the taylor-series expansion
    intKrnl.power = i_timeStepWidth;

    intKrnl.execute0();

    for (unsigned der = 1; der < Config::ConvergenceOrder; ++der) {
      krnl.execute(der);

      // update scalar for this derivative
      intKrnl.power *= i_timeStepWidth / RealT(der + 1);
      intKrnl.execute(der);
    }

    // TODO(Lukas) Implement!
    // Compute integrated displacement over time step if needed.
  }

  void flopsAder(unsigned int& o_nonZeroFlops, unsigned int& o_hardwareFlops) {
    // reset flops
    o_nonZeroFlops = 0;
    o_hardwareFlops = 0;

    // initialization
    o_nonZeroFlops += Yateto<Config>::Kernel::derivativeTaylorExpansion::nonZeroFlops(0);
    o_hardwareFlops += Yateto<Config>::Kernel::derivativeTaylorExpansion::hardwareFlops(0);

    // interate over derivatives
    for (unsigned der = 1; der < Config::ConvergenceOrder; ++der) {
      o_nonZeroFlops += Yateto<Config>::Kernel::derivative::nonZeroFlops(der);
      o_hardwareFlops += Yateto<Config>::Kernel::derivative::hardwareFlops(der);

      // update of time integrated DOFs
      o_nonZeroFlops += Yateto<Config>::Kernel::derivativeTaylorExpansion::nonZeroFlops(der);
      o_hardwareFlops += Yateto<Config>::Kernel::derivativeTaylorExpansion::hardwareFlops(der);
    }
  }

  unsigned bytesAder() {
    unsigned reals = 0;

    // DOFs load, tDOFs load, tDOFs write
    reals += Yateto<Config>::Tensor::Q::size() + Yateto<Config>::Tensor::Qane::size() +
             2 * Yateto<Config>::Tensor::I::size() + 2 * Yateto<Config>::Tensor::Iane::size();
    // star matrices, source matrix
    reals += yateto::computeFamilySize<typename Yateto<Config>::Tensor::star>() +
             Yateto<Config>::Tensor::w::size() + Yateto<Config>::Tensor::W::size() +
             Yateto<Config>::Tensor::E::size();

    /// \todo incorporate derivatives

    return reals * sizeof(RealT);
  }

  void computeIntegral(double i_expansionPoint,
                       double i_integrationStart,
                       double i_integrationEnd,
                       RealT const* i_timeDerivatives,
                       RealT o_timeIntegrated[ConfigConstants<Config>::TensorSizeI]) {

    /*
     * assert alignments.
     */
    assert(((uintptr_t)i_timeDerivatives) % Alignment == 0);
    assert(((uintptr_t)o_timeIntegrated) % Alignment == 0);

    // compute lengths of integration intervals
    RealT l_deltaTLower = i_integrationStart - i_expansionPoint;
    RealT l_deltaTUpper = i_integrationEnd - i_expansionPoint;

    // initialization of scalars in the taylor series expansion (0th term)
    RealT l_firstTerm = static_cast<RealT>(1.0);
    RealT l_secondTerm = static_cast<RealT>(1.0);
    RealT l_factorial = static_cast<RealT>(1.0);

    typename Yateto<Config>::Kernel::derivativeTaylorExpansionEla intKrnl;
    intKrnl.I = o_timeIntegrated;
    RealT const* der = i_timeDerivatives;
    for (unsigned i = 0; i < yateto::numFamilyMembers<Yateto<Config>::Tensor::dQ>(); ++i) {
      intKrnl.dQ(i) = der;
      der += Yateto<Config>::Tensor::dQ::size(i);
    }

    // iterate over time derivatives
    for (int der = 0; der < Config::ConvergenceOrder; ++der) {
      l_firstTerm *= l_deltaTUpper;
      l_secondTerm *= l_deltaTLower;
      l_factorial *= static_cast<RealT>(der + 1);

      intKrnl.power = l_firstTerm - l_secondTerm;
      intKrnl.power /= l_factorial;

      intKrnl.execute(der);
    }
  }

  void computeTaylorExpansion(RealT time,
                              RealT expansionPoint,
                              RealT const* timeDerivatives,
                              RealT timeEvaluated[Yateto<Config>::Tensor::Q::size()]) {
    /*
     * assert alignments.
     */
    assert(((uintptr_t)timeDerivatives) % Alignment == 0);
    assert(((uintptr_t)timeEvaluated) % Alignment == 0);

    // assert that this is a forward evaluation in time
    assert(time >= expansionPoint);

    RealT deltaT = time - expansionPoint;

    static_assert(Yateto<Config>::Tensor::I::size() == Yateto<Config>::Tensor::Q::size(),
                  "Sizes of tensors I and Q must match");

    typename Yateto<Config>::Kernel::derivativeTaylorExpansionEla intKrnl;
    intKrnl.I = timeEvaluated;
    RealT const* der = timeDerivatives;
    for (unsigned i = 0; i < yateto::numFamilyMembers<Yateto<Config>::Tensor::dQ>(); ++i) {
      intKrnl.dQ(i) = der;
      der += Yateto<Config>::Tensor::dQ::size(i);
    }
    intKrnl.power = 1.0;

    // iterate over time derivatives
    for (int derivative = 0; derivative < Config::ConvergenceOrder; ++derivative) {
      intKrnl.execute(derivative);
      intKrnl.power *= deltaT / RealT(derivative + 1);
    }
  }

  void flopsTaylorExpansion(long long& nonZeroFlops, long long& hardwareFlops) {
    // reset flops
    nonZeroFlops = 0;
    hardwareFlops = 0;

    // interate over derivatives
    for (unsigned der = 0; der < Config::ConvergenceOrder; ++der) {
      nonZeroFlops += Yateto<Config>::Kernel::derivativeTaylorExpansionEla::nonZeroFlops(der);
      hardwareFlops += Yateto<Config>::Kernel::derivativeTaylorExpansionEla::hardwareFlops(der);
    }
  }
};
} // namespace seissol::waveprop::kernel::time

#endif
