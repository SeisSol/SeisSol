/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Stephanie Wollherr (wollherr AT geophysik.uni-muenchen.de, https://www.geophysik.uni-muenchen.de/Members/wollherr)
 *
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
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
 * Plasticity kernel of SeisSol.
 **/

#include "Plasticity.h"

#include <cstring>
#include <algorithm>
#include <cmath>
#include <generated_code/kernel.h>
#include <generated_code/init.h>
#include "common.hpp"

#ifdef ACL_DEVICE
#include "device.h"
#include "DeviceAux/PlasticityAux.h"
using namespace device;
#endif

namespace seissol::kernels {
  unsigned Plasticity::computePlasticity(double oneMinusIntegratingFactor,
                                         double timeStepWidth,
                                         double T_v,
                                         GlobalData const *global,
                                         PlasticityData const *plasticityData,
                                         real degreesOfFreedom[tensor::Q::size()],
                                         real *pstrain) {
    assert(reinterpret_cast<uintptr_t>(degreesOfFreedom) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(global->vandermondeMatrix) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(global->vandermondeMatrixInverse) % ALIGNMENT == 0);

    real QStressNodal[tensor::QStressNodal::size()] __attribute__((aligned(ALIGNMENT)));
    real meanStress[tensor::meanStress::size()] __attribute__((aligned(ALIGNMENT)));
    real secondInvariant[tensor::secondInvariant::size()] __attribute__((aligned(ALIGNMENT)));
    real tau[tensor::secondInvariant::size()] __attribute__((aligned(ALIGNMENT)));
    real taulim[tensor::meanStress::size()] __attribute__((aligned(ALIGNMENT)));
    real yieldFactor[tensor::yieldFactor::size()] __attribute__((aligned(ALIGNMENT)));
    real dudt_pstrain[7];

    static_assert(tensor::secondInvariant::size() == tensor::meanStress::size(),
                  "Second invariant tensor and mean stress tensor must be of the same size().");
    static_assert(tensor::yieldFactor::size() <= tensor::meanStress::size(),
                  "Yield factor tensor must be smaller than mean stress tensor.");

    //copy dofs for later comparison, only first dof of stresses required
    // @todo multiple sims
    real prev_degreesOfFreedom[6];
    for (unsigned q = 0; q < 6; ++q) {
      prev_degreesOfFreedom[q] = degreesOfFreedom[q * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    }

    /* Convert modal to nodal and add sigma0.
     * Stores s_{ij} := sigma_{ij} + sigma0_{ij} for every node.
     * sigma0 is constant */
    kernel::plConvertToNodal m2nKrnl;
    m2nKrnl.v = global->vandermondeMatrix;
    m2nKrnl.QStress = degreesOfFreedom;
    m2nKrnl.QStressNodal = QStressNodal;
    m2nKrnl.replicateInitialLoading = init::replicateInitialLoading::Values;
    m2nKrnl.initialLoading = plasticityData->initialLoading;
    m2nKrnl.execute();

    // Computes m = s_{ii} / 3.0 for every node
    kernel::plComputeMean cmKrnl;
    cmKrnl.meanStress = meanStress;
    cmKrnl.QStressNodal = QStressNodal;
    cmKrnl.selectBulkAverage = init::selectBulkAverage::Values;
    cmKrnl.execute();

    /* Compute s_{ij} := s_{ij} - m delta_{ij},
     * where delta_{ij} = 1 if i == j else 0.
     * Thus, s_{ij} contains the deviatoric stresses. */
    kernel::plSubtractMean smKrnl;
    smKrnl.meanStress = meanStress;
    smKrnl.QStressNodal = QStressNodal;
    smKrnl.selectBulkNegative = init::selectBulkNegative::Values;
    smKrnl.execute();

    // Compute I_2 = 0.5 s_{ij} s_ji for every node
    kernel::plComputeSecondInvariant siKrnl;
    siKrnl.secondInvariant = secondInvariant;
    siKrnl.QStressNodal = QStressNodal;
    siKrnl.weightSecondInvariant = init::weightSecondInvariant::Values;
    siKrnl.execute();

    // tau := sqrt(I_2) for every node
    for (unsigned ip = 0; ip < tensor::secondInvariant::size(); ++ip) {
      tau[ip] = sqrt(secondInvariant[ip]);
    }

    // Compute tau_c for every node
    for (unsigned ip = 0; ip < tensor::meanStress::size(); ++ip) {
      taulim[ip] = std::max((real) 0.0, plasticityData->cohesionTimesCosAngularFriction -
                                        meanStress[ip] * plasticityData->sinAngularFriction);
    }

    bool adjust = false;
    for (unsigned ip = 0; ip < tensor::yieldFactor::size(); ++ip) {
      // Compute yield := (t_c / tau - 1) r for every node,
      // where r = 1 - exp(-timeStepWidth / T_v)
      if (tau[ip] > taulim[ip]) {
        adjust = true;
        yieldFactor[ip] = (taulim[ip] / tau[ip] - 1.0) * oneMinusIntegratingFactor;
      } else {
        yieldFactor[ip] = 0.0;
      }
    }

    if (adjust) {
      /**
       * Compute sigma_{ij} := sigma_{ij} + yield s_{ij} for every node
       * and store as modal basis.
       *
       * Remark: According to Wollherr et al., the update formula (13) should be
       *
       * sigmaNew_{ij} := f^* s_{ij} + m delta_{ij} - sigma0_{ij}
       *
       * where f^* = r tau_c / tau + (1 - r) = 1 + yield. Adding 0 to (13) gives
       *
       * sigmaNew_{ij} := f^* s_{ij} + m delta_{ij} - sigma0_{ij}
       *                  + sigma_{ij} + sigma0_{ij} - sigma_{ij} - sigma0_{ij}
       *                = f^* s_{ij} + sigma_{ij} - s_{ij}
       *                = sigma_{ij} + (f^* - 1) s_{ij}
       *                = sigma_{ij} + yield s_{ij}
       */
      kernel::plAdjustStresses adjKrnl;
      adjKrnl.QStress = degreesOfFreedom;
      adjKrnl.vInv = global->vandermondeMatrixInverse;
      adjKrnl.QStressNodal = QStressNodal;
      adjKrnl.yieldFactor = yieldFactor;
      adjKrnl.execute();

      // calculate plastic strain with first dof only (for now)
      for (unsigned q = 0; q < 6; ++q) {
        /**
         * Equation (10) from Wollherr et al.:
         *
         * d/dt strain_{ij} = (sigma_{ij} + sigma0_{ij} - P_{ij}(sigma)) / (2mu T_v)
         *
         * where (11)
         *
         * P_{ij}(sigma) = { tau_c/tau s_{ij} + m delta_{ij}         if     tau >= taulim
         *                 { sigma_{ij} + sigma0_{ij}                else
         *
         * Thus,
         *
         * d/dt strain_{ij} = { (1 - tau_c/tau) / (2mu T_v) s_{ij}   if     tau >= taulim
         *                    { 0                                    else
         *
         * Consider tau >= taulim first. We have (1 - tau_c/tau) = -yield / r. Therefore,
         *
         * d/dt strain_{ij} = -1 / (2mu T_v r) yield s_{ij}
         *                  = -1 / (2mu T_v r) (sigmaNew_{ij} - sigma_{ij})
         *                  = (sigma_{ij} - sigmaNew_{ij}) / (2mu T_v r)
         *
         * If tau < taulim, then sigma_{ij} - sigmaNew_{ij} = 0.
         */
        real factor = plasticityData->mufactor / (T_v * oneMinusIntegratingFactor);
        dudt_pstrain[q] = factor *
            (prev_degreesOfFreedom[q] - degreesOfFreedom[q * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS]);
        // Integrate with explicit Euler
        pstrain[q] += timeStepWidth * dudt_pstrain[q];
      }

      // eta := int_0^t sqrt(0.5 dstrain_{ij}/dt dstrain_{ij}/dt) dt
      // Approximate with eta += timeStepWidth * sqrt(0.5 dstrain_{ij}/dt dstrain_{ij}/dt)
      pstrain[6] += timeStepWidth * sqrt(0.5 * (dudt_pstrain[0] * dudt_pstrain[0] + dudt_pstrain[1] * dudt_pstrain[1]
                                                + dudt_pstrain[2] * dudt_pstrain[2]) + dudt_pstrain[3] * dudt_pstrain[3]
                                         + dudt_pstrain[4] * dudt_pstrain[4] + dudt_pstrain[5] * dudt_pstrain[5]);

      return 1;
    }

    return 0;
  }

  unsigned Plasticity::computePlasticityBatched(double oneMinusIntegratingFactor,
                                                double timeStepWidth,
                                                double T_v,
                                                GlobalData const *global,
                                                initializers::recording::ConditionalBatchTableT &table,
                                                PlasticityData *plasticity) {
#ifdef ACL_DEVICE
    static_assert(tensor::Q::Shape[0] == tensor::QStressNodal::Shape[0],
                  "modal and nodal dofs must have the same leading dimensions");
    static_assert(tensor::Q::Shape[0] == tensor::v::Shape[0],
                  "modal dofs and vandermonde matrix must hage the same leading dimensions");

    DeviceInstance &device = DeviceInstance::getInstance();
    ConditionalKey key(*KernelNames::Plasticity);

    if (table.find(key) != table.end()) {
      unsigned stackMemCounter{0};
      BatchTable& entry = table[key];
      const size_t numElements = (entry.content[*EntityId::Dofs])->getSize();
      real** modalStressTensors = (entry.content[*EntityId::Dofs])->getPointers();
      real** nodalStressTensors = (entry.content[*EntityId::NodalStressTensor])->getPointers();

      const size_t firsModesSize = NUM_STREESS_COMPONENTS * numElements * sizeof(real);
      real *firsModes = reinterpret_cast<real*>(device.api->getStackMemory(firsModesSize));
      ++stackMemCounter;

      device::aux::plasticity::saveFirstModes(firsModes,
                                              const_cast<const real**>(modalStressTensors),
                                              numElements);

      assert(global->replicateStresses != nullptr && "replicateStresses has not been initialized");
      real** initLoad = (entry.content[*EntityId::InitialLoad])->getPointers();
      kernel::gpu_plConvertToNodal m2nKrnl;
      m2nKrnl.v = global->vandermondeMatrix;
      m2nKrnl.QStress = const_cast<const real**>(modalStressTensors);
      m2nKrnl.QStressNodal = nodalStressTensors;
      m2nKrnl.replicateInitialLoadingM = global->replicateStresses;
      m2nKrnl.initialLoadingM = const_cast<const real**>(initLoad);
      m2nKrnl.numElements = numElements;

      const size_t MAX_TMP_MEM = m2nKrnl.TmpMaxMemRequiredInBytes * numElements;
      real *tmpMem = (real*)(device.api->getStackMemory(MAX_TMP_MEM));
      ++stackMemCounter;

      m2nKrnl.linearAllocator.initialize(tmpMem);
      m2nKrnl.execute();
      m2nKrnl.linearAllocator.free();

      int *isAdjustableVector =
          reinterpret_cast<int*>(device.api->getStackMemory(numElements * sizeof(int)));
      ++stackMemCounter;

      device::aux::plasticity::adjustDeviatoricTensors(nodalStressTensors,
                                                       isAdjustableVector,
                                                       plasticity,
                                                       oneMinusIntegratingFactor,
                                                       numElements);

      unsigned numAdjustedDofs = device.algorithms.reduceVector(isAdjustableVector,
                                                                numElements,
                                                                ::device::ReductionType::Add);

      // apply stress adjustment
      device::aux::plasticity::adjustModalStresses(modalStressTensors,
                                                   const_cast<const real **>(nodalStressTensors),
                                                   global->vandermondeMatrixInverse,
                                                   isAdjustableVector,
                                                   numElements);

      // compute Pstrains
      real **pstrains = entry.content[*EntityId::Pstrains]->getPointers();
      device::aux::plasticity::computePstrains(pstrains,
                                               isAdjustableVector,
                                               const_cast<const real **>(modalStressTensors),
                                               firsModes,
                                               plasticity,
                                               oneMinusIntegratingFactor,
                                               timeStepWidth,
                                               T_v,
                                               numElements);

      // NOTE: Temp memory must be properly clean after using negative signed integers
      // This kind of memory is mainly used for floating-point numbers. Negative signed ints might corrupt
      // the most significant bits. We came to this conclusion by our first-hand experience
      device.algorithms.fillArray(reinterpret_cast<char*>(isAdjustableVector),
                                  static_cast<char>(0),
                                  numElements * sizeof(int));

      for (unsigned i = 0; i < stackMemCounter; ++i) {
        device.api->popStackMemory();
      }
      return numAdjustedDofs;
    }
    else {
      return 0;
    }

#else
    assert(false && "no implementation provided");
    return 0;
#endif // ACL_DEVICE
  }

  void Plasticity::flopsPlasticity(long long &o_NonZeroFlopsCheck,
                                   long long &o_HardwareFlopsCheck,
                                   long long &o_NonZeroFlopsYield,
                                   long long &o_HardwareFlopsYield) {
    // reset flops
    o_NonZeroFlopsCheck = 0;
    o_HardwareFlopsCheck = 0;
    o_NonZeroFlopsYield = 0;
    o_HardwareFlopsYield = 0;

    // flops from checking, i.e. outside if (adjust) {}
    o_NonZeroFlopsCheck += kernel::plConvertToNodal::NonZeroFlops;
    o_HardwareFlopsCheck += kernel::plConvertToNodal::HardwareFlops;

    // compute mean stress
    o_NonZeroFlopsCheck += kernel::plComputeMean::NonZeroFlops;
    o_HardwareFlopsCheck += kernel::plComputeMean::HardwareFlops;

    // subtract mean stress
    o_NonZeroFlopsCheck += kernel::plSubtractMean::NonZeroFlops;
    o_HardwareFlopsCheck += kernel::plSubtractMean::HardwareFlops;

    // compute second invariant
    o_NonZeroFlopsCheck += kernel::plComputeSecondInvariant::NonZeroFlops;
    o_HardwareFlopsCheck += kernel::plComputeSecondInvariant::HardwareFlops;

    // compute taulim (1 add, 1 mul, max NOT counted)
    o_NonZeroFlopsCheck += 2 * tensor::meanStress::size();
    o_HardwareFlopsCheck += 2 * tensor::meanStress::size();

    // check for yield (NOT counted, as it would require counting the number of yielding points)

    // flops from plastic yielding, i.e. inside if (adjust) {}
    o_NonZeroFlopsYield += kernel::plAdjustStresses::NonZeroFlops;
    o_HardwareFlopsYield += kernel::plAdjustStresses::HardwareFlops;
  }
} // namespace seissol::kernels
