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

#ifdef _OPENMP
#  include <omp.h>
#endif

void seissol::initializers::initializeGlobalData(GlobalData& globalData, memory::ManagedAllocator& memoryAllocator, enum seissol::memory::Memkind memkind)
{
  // DG global matrices
  real* globalMatrixMem = static_cast<real*>(memoryAllocator.allocateMemory( seissol::model::globalMatrixOffsets[seissol::model::numGlobalMatrices] * sizeof(real), PAGESIZE_HEAP, memkind ));
  for (unsigned matrix = 0; matrix < seissol::model::numGlobalMatrices; ++matrix) {
    memcpy(
      &globalMatrixMem[ seissol::model::globalMatrixOffsets[matrix] ],
      seissol::model::globalMatrixValues[matrix],
      (seissol::model::globalMatrixOffsets[matrix+1] - seissol::model::globalMatrixOffsets[matrix]) * sizeof(real)
    );
  }
  for (unsigned transposedStiffness = 0; transposedStiffness < 3; ++transposedStiffness) {
    globalData.stiffnessMatricesTransposed[transposedStiffness] = &globalMatrixMem[ seissol::model::globalMatrixOffsets[transposedStiffness] ];
  }
  for (unsigned stiffness = 0; stiffness < 3; ++stiffness) {
    globalData.stiffnessMatrices[stiffness] = &globalMatrixMem[ seissol::model::globalMatrixOffsets[3 + stiffness] ];
  }
  for (unsigned cob = 0; cob < 4; ++cob) {
    globalData.changeOfBasisMatrices[cob] = &globalMatrixMem[ seissol::model::globalMatrixOffsets[6 + cob] ];
  }
  for (unsigned cob = 0; cob < 4; ++cob) {
    globalData.neighbourChangeOfBasisMatricesTransposed[cob] = &globalMatrixMem[ seissol::model::globalMatrixOffsets[10 + cob] ];
  }
  for (unsigned cob = 0; cob < 4; ++cob) {
    globalData.localChangeOfBasisMatricesTransposed[cob] = &globalMatrixMem[ seissol::model::globalMatrixOffsets[14 + cob] ];
  }
  for (unsigned flux = 0; flux < 3; ++flux) {
    globalData.neighbourFluxMatrices[flux] = &globalMatrixMem[ seissol::model::globalMatrixOffsets[18 + flux] ];
  }

  // @TODO Integrate this step into the code generator
  for (unsigned transposedStiffness = 0; transposedStiffness < 3; ++transposedStiffness) {
    real* matrix = &globalMatrixMem[ seissol::model::globalMatrixOffsets[transposedStiffness] ];
    for (unsigned i = 0; i < seissol::model::globalMatrixOffsets[transposedStiffness+1]-seissol::model::globalMatrixOffsets[transposedStiffness]; ++i) {
      matrix[i] *= -1.0;
    }
  }

  // Dynamic Rupture global matrices
  real* drGlobalMatrixMem = static_cast<real*>(memoryAllocator.allocateMemory( seissol::model::dr_globalMatrixOffsets[seissol::model::dr_numGlobalMatrices] * sizeof(real), PAGESIZE_HEAP, memkind ));
  for (unsigned matrix = 0; matrix < seissol::model::dr_numGlobalMatrices; ++matrix) {
    memcpy(
      &drGlobalMatrixMem[ seissol::model::dr_globalMatrixOffsets[matrix] ],
      seissol::model::dr_globalMatrixValues[matrix],
      (seissol::model::dr_globalMatrixOffsets[matrix+1] - seissol::model::dr_globalMatrixOffsets[matrix]) * sizeof(real)
    );
  }
  for (unsigned face = 0; face < 4; ++face) {
    for (unsigned h = 0; h < 4; ++h) {
      globalData.nodalFluxMatrices[face][h] = &drGlobalMatrixMem[ seissol::model::dr_globalMatrixOffsets[4*face+h] ];
      globalData.faceToNodalMatrices[face][h] = &drGlobalMatrixMem[ seissol::model::dr_globalMatrixOffsets[16 + 4*face+h] ];
    }
  }

  real* plasticityGlobalMatrixMem = static_cast<real*>(memoryAllocator.allocateMemory( seissol::model::plasticity_globalMatrixOffsets[seissol::model::plasticity_numGlobalMatrices] * sizeof(real), PAGESIZE_HEAP, memkind ));
  for (unsigned matrix = 0; matrix < seissol::model::plasticity_numGlobalMatrices; ++matrix) {
    memcpy(
      &plasticityGlobalMatrixMem[ seissol::model::plasticity_globalMatrixOffsets[matrix] ],
      seissol::model::plasticity_globalMatrixValues[matrix],
      (seissol::model::plasticity_globalMatrixOffsets[matrix+1] - seissol::model::plasticity_globalMatrixOffsets[matrix]) * sizeof(real)
    );
  }
  globalData.vandermondeMatrix = &plasticityGlobalMatrixMem[ seissol::model::plasticity_globalMatrixOffsets[0] ];
  globalData.vandermondeMatrixInverse = &plasticityGlobalMatrixMem[ seissol::model::plasticity_globalMatrixOffsets[1] ];
  
  // thread-local LTS integration buffers  
  int l_numberOfThreads = 1;
#ifdef _OPENMP
  l_numberOfThreads = omp_get_max_threads();
#endif
  real* integrationBufferLTS = (real*) memoryAllocator.allocateMemory( l_numberOfThreads*(4*NUMBER_OF_ALIGNED_DOFS)*sizeof(real), PAGESIZE_STACK, memkind ) ;

  // initialize w.r.t. NUMA
#ifdef _OPENMP
  #pragma omp parallel
  {
    size_t l_threadOffset = omp_get_thread_num()*(4*NUMBER_OF_ALIGNED_DOFS);
#else
    size_t l_threadOffset = 0;
#endif
    for ( unsigned int l_dof = 0; l_dof < (4*NUMBER_OF_ALIGNED_DOFS); l_dof++ ) {
      integrationBufferLTS[l_dof + l_threadOffset] = (real)0.0;
    }
#ifdef _OPENMP
  }
#endif
  
  globalData.integrationBufferLTS = integrationBufferLTS;
}
