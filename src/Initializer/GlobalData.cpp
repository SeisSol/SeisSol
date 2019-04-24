/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 **/

#include "GlobalData.h"
#include <generated_code/init.h>
#include <yateto.h>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace init = seissol::init;

void seissol::initializers::initializeGlobalData(GlobalData& globalData, memory::ManagedAllocator& memoryAllocator, enum seissol::memory::Memkind memkind)
{
  // We ensure that global matrices always start at an aligned memory address,
  // such that mixed cases with aligned and non-aligned global matrices do also work.  

  unsigned globalMatrixMemSize = 0;
  globalMatrixMemSize += yateto::computeFamilySize<init::kDivM>(yateto::alignedReals<real>(ALIGNMENT));
  globalMatrixMemSize += yateto::computeFamilySize<init::kDivMT>(yateto::alignedReals<real>(ALIGNMENT));
  globalMatrixMemSize += yateto::computeFamilySize<init::rDivM>(yateto::alignedReals<real>(ALIGNMENT));
  globalMatrixMemSize += yateto::computeFamilySize<init::rT>(yateto::alignedReals<real>(ALIGNMENT));
  globalMatrixMemSize += yateto::computeFamilySize<init::fMrT>(yateto::alignedReals<real>(ALIGNMENT));
  globalMatrixMemSize += yateto::computeFamilySize<init::fP>(yateto::alignedReals<real>(ALIGNMENT));
  globalMatrixMemSize += yateto::alignedUpper(tensor::evalAtQP::size(),  yateto::alignedReals<real>(ALIGNMENT));
  globalMatrixMemSize += yateto::alignedUpper(tensor::projectQP::size(), yateto::alignedReals<real>(ALIGNMENT));
  
  real* globalMatrixMem = static_cast<real*>(memoryAllocator.allocateMemory( globalMatrixMemSize * sizeof(real), PAGESIZE_HEAP, memkind ));

  real* globalMatrixMemPtr = globalMatrixMem;
  yateto::copyFamilyToMemAndSetPtr<init::kDivMT, real>(globalMatrixMemPtr, globalData.stiffnessMatricesTransposed, ALIGNMENT);
  yateto::copyFamilyToMemAndSetPtr<init::kDivM,  real>(globalMatrixMemPtr, globalData.stiffnessMatrices, ALIGNMENT);
  yateto::copyFamilyToMemAndSetPtr<init::rDivM,  real>(globalMatrixMemPtr, globalData.changeOfBasisMatrices, ALIGNMENT);
  yateto::copyFamilyToMemAndSetPtr<init::rT,     real>(globalMatrixMemPtr, globalData.neighbourChangeOfBasisMatricesTransposed, ALIGNMENT);
  yateto::copyFamilyToMemAndSetPtr<init::fMrT,   real>(globalMatrixMemPtr, globalData.localChangeOfBasisMatricesTransposed, ALIGNMENT);
  yateto::copyFamilyToMemAndSetPtr<init::fP,     real>(globalMatrixMemPtr, globalData.neighbourFluxMatrices, ALIGNMENT);
  yateto::copyTensorToMemAndSetPtr<init::evalAtQP,     real>(globalMatrixMemPtr, globalData.evalAtQPMatrix, ALIGNMENT);
  yateto::copyTensorToMemAndSetPtr<init::projectQP,    real>(globalMatrixMemPtr, globalData.projectQPMatrix, ALIGNMENT);
  
  assert(globalMatrixMemPtr == globalMatrixMem + globalMatrixMemSize);

  // @TODO Integrate this step into the code generator
  for (unsigned transposedStiffness = 0; transposedStiffness < 3; ++transposedStiffness) {
    real* matrix = const_cast<real*>(globalData.stiffnessMatricesTransposed(transposedStiffness));
    for (unsigned i = 0; i < init::kDivMT::size(transposedStiffness); ++i) {
      matrix[i] *= -1.0;
    }
  }

  // Dynamic Rupture global matrices
  unsigned drGlobalMatrixMemSize = 0;
  drGlobalMatrixMemSize += yateto::computeFamilySize<init::V3mTo2nTWDivM>(yateto::alignedReals<real>(ALIGNMENT));
  drGlobalMatrixMemSize += yateto::computeFamilySize<init::V3mTo2n>(yateto::alignedReals<real>(ALIGNMENT));
  
  real* drGlobalMatrixMem = static_cast<real*>(memoryAllocator.allocateMemory( drGlobalMatrixMemSize  * sizeof(real), PAGESIZE_HEAP, memkind ));
  
  real* drGlobalMatrixMemPtr = drGlobalMatrixMem;
  yateto::copyFamilyToMemAndSetPtr<init::V3mTo2nTWDivM, real>(drGlobalMatrixMemPtr, globalData.nodalFluxMatrices, ALIGNMENT);
  yateto::copyFamilyToMemAndSetPtr<init::V3mTo2n,       real>(drGlobalMatrixMemPtr, globalData.faceToNodalMatrices, ALIGNMENT);
  
  assert(drGlobalMatrixMemPtr == drGlobalMatrixMem + drGlobalMatrixMemSize);

  // Plasticity global matrices
  unsigned plasticityGlobalMatrixMemSize = 0;
  plasticityGlobalMatrixMemSize += yateto::alignedUpper(tensor::v::size(),    yateto::alignedReals<real>(ALIGNMENT));
  plasticityGlobalMatrixMemSize += yateto::alignedUpper(tensor::vInv::size(), yateto::alignedReals<real>(ALIGNMENT));

  real* plasticityGlobalMatrixMem = static_cast<real*>(memoryAllocator.allocateMemory( plasticityGlobalMatrixMemSize * sizeof(real), PAGESIZE_HEAP, memkind ));
  
  real* plasticityGlobalMatrixMemPtr = plasticityGlobalMatrixMem;
  yateto::copyTensorToMemAndSetPtr<init::v,    real>(plasticityGlobalMatrixMemPtr, globalData.vandermondeMatrix, ALIGNMENT);
  yateto::copyTensorToMemAndSetPtr<init::vInv, real>(plasticityGlobalMatrixMemPtr, globalData.vandermondeMatrixInverse, ALIGNMENT);
  
  assert(plasticityGlobalMatrixMemPtr == plasticityGlobalMatrixMem + plasticityGlobalMatrixMemSize);
  
  // thread-local LTS integration buffers  
  int l_numberOfThreads = 1;
#ifdef _OPENMP
  l_numberOfThreads = omp_get_max_threads();
#endif
  real* integrationBufferLTS = (real*) memoryAllocator.allocateMemory( l_numberOfThreads*(4*tensor::I::size())*sizeof(real), PAGESIZE_STACK, memkind ) ;

  // initialize w.r.t. NUMA
#ifdef _OPENMP
  #pragma omp parallel
  {
    size_t l_threadOffset = omp_get_thread_num()*(4*tensor::I::size());
#else
    size_t l_threadOffset = 0;
#endif
    for ( unsigned int l_dof = 0; l_dof < (4*tensor::I::size()); l_dof++ ) {
      integrationBufferLTS[l_dof + l_threadOffset] = (real)0.0;
    }
#ifdef _OPENMP
  }
#endif
  
  globalData.integrationBufferLTS = integrationBufferLTS;
}
