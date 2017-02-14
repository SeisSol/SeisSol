/*
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


void computeAderIntegration() {
#ifdef _OPENMP
  #pragma omp parallel 
  {
#if NUMBER_OF_GLOBAL_DATA_COPIES>1
  //GlobalData* l_globalData = m_globalDataArray[(omp_get_thread_num()/NUMBER_OF_COMPACT_THREADS_PER_GLOBAL_DATA_COPY)%NUMBER_OF_GLOBAL_DATA_COPIES];
  GlobalData* l_globalData = m_globalDataArray[0];
#else
  GlobalData* l_globalData = m_globalData;
#endif
  #pragma omp for schedule(static)
#else
  GlobalData* l_globalData = m_globalData;
#endif
  for( unsigned int l_cell = 0; l_cell < m_cells->numberOfCells; l_cell++ ) {
    m_timeKernel.computeAder(              m_timeStepWidthSimulation,
                                           l_globalData,
                                           &m_cellData->localIntegration[l_cell],
                                           m_cells->dofs[l_cell],
                                           m_cells->buffers[l_cell],
                                           m_cells->derivatives[l_cell] );
  }
#ifdef _OPENMP
  }
#endif
}

void computeLocalWithoutAderIntegration() {
#ifdef _OPENMP
  #pragma omp parallel 
  {
#if NUMBER_OF_GLOBAL_DATA_COPIES>1
  //GlobalData* l_globalData = m_globalDataArray[(omp_get_thread_num()/NUMBER_OF_COMPACT_THREADS_PER_GLOBAL_DATA_COPY)%NUMBER_OF_GLOBAL_DATA_COPIES];
  GlobalData* l_globalData = m_globalDataArray[0];
#else
  GlobalData* l_globalData = m_globalData;
#endif
  #pragma omp for schedule(static)
#else
  GlobalData* l_globalData = m_globalData;
#endif
  for( unsigned int l_cell = 0; l_cell < m_cells->numberOfCells; l_cell++ ) {
    m_localKernel.computeIntegral(  m_cellInformation[l_cell].faceTypes,
                                    l_globalData,
                                    &m_cellData->localIntegration[l_cell],
                                    m_cells->buffers[l_cell],
                                    m_cells->dofs[l_cell] );
  }
#ifdef _OPENMP
  }
#endif
}

void computeLocalIntegration() {
#ifdef _OPENMP
  #pragma omp parallel
  {
#if NUMBER_OF_GLOBAL_DATA_COPIES>1
  //GlobalData* l_globalData = m_globalDataArray[(omp_get_thread_num()/NUMBER_OF_COMPACT_THREADS_PER_GLOBAL_DATA_COPY)%NUMBER_OF_GLOBAL_DATA_COPIES];
  GlobalData* l_globalData = m_globalDataArray[0];
#else
  GlobalData* l_globalData = m_globalData;
#endif
  #pragma omp for schedule(static)
#else
  GlobalData* l_globalData = m_globalData;
#endif
  for( unsigned int l_cell = 0; l_cell < m_cells->numberOfCells; l_cell++ ) {
    m_timeKernel.computeAder(      (double)m_timeStepWidthSimulation,
                                           l_globalData,
                                           &m_cellData->localIntegration[l_cell],
                                           m_cells->dofs[l_cell],
                                           m_cells->buffers[l_cell],
                                           m_cells->derivatives[l_cell] );

    m_localKernel.computeIntegral(        m_cellInformation[l_cell].faceTypes,
                                          l_globalData,
                                           &m_cellData->localIntegration[l_cell],
                                           m_cells->buffers[l_cell],
                                           m_cells->dofs[l_cell] );
  }
#ifdef _OPENMP
  }
#endif
}

void computeNeighboringIntegration() {
  real  l_integrationBuffer[4][NUMBER_OF_ALIGNED_DOFS] __attribute__((aligned(4096)));
  real *l_timeIntegrated[4];
#ifdef ENABLE_MATRIX_PREFETCH
  real *l_faceNeighbors_prefetch[4];
#endif

#ifdef _OPENMP
#  ifdef ENABLE_MATRIX_PREFETCH
  #pragma omp parallel private(l_integrationBuffer, l_timeIntegrated, l_faceNeighbors_prefetch)
#  else
  #pragma omp parallel private(l_integrationBuffer, l_timeIntegrated)
#  endif
  {
#if NUMBER_OF_THREADS_PER_GLOBALDATA_COPY < 512
  GlobalData* l_globalData = m_globalDataArray[omp_get_thread_num()/NUMBER_OF_THREADS_PER_GLOBALDATA_COPY];
#else
  GlobalData* l_globalData = m_globalData;
#endif
  #pragma omp for schedule(static)
#else
  GlobalData* l_globalData = m_globalData;
#endif
  for( int l_cell = 0; l_cell < m_cells->numberOfCells; l_cell++ ) {
    seissol::kernels::TimeCommon::computeIntegrals(m_timeKernel,
                                              m_cellInformation[l_cell].ltsSetup,
                                               m_cellInformation[l_cell].faceTypes,
                                               0.0,
                                       (double)m_timeStepWidthSimulation,
                                               m_cells->faceNeighbors[l_cell],
                                               l_integrationBuffer,
                                               l_timeIntegrated );

#ifdef ENABLE_MATRIX_PREFETCH
#pragma message("the current prefetch structure (flux matrices and tDOFs is tuned for higher order and shouldn't be harmful for lower orders")
    l_faceNeighbors_prefetch[0] = (m_cellInformation[l_cell].faceTypes[1] != dynamicRupture) ? m_cells->faceNeighbors[l_cell][1] : m_cells->drMapping[l_cell][1].godunov;
    l_faceNeighbors_prefetch[1] = (m_cellInformation[l_cell].faceTypes[2] != dynamicRupture) ? m_cells->faceNeighbors[l_cell][2] : m_cells->drMapping[l_cell][2].godunov;
    l_faceNeighbors_prefetch[2] = (m_cellInformation[l_cell].faceTypes[3] != dynamicRupture) ? m_cells->faceNeighbors[l_cell][3] : m_cells->drMapping[l_cell][3].godunov;

    // fourth face's prefetches
    if (l_cell < (m_cells->numberOfCells-1) ) {
      l_faceNeighbors_prefetch[3] = (m_cellInformation[l_cell+1].faceTypes[0] != dynamicRupture) ? m_cells->faceNeighbors[l_cell+1][0] : m_cells->drMapping[l_cell+1][0].godunov;
    } else {
      l_faceNeighbors_prefetch[3] = m_cells->faceNeighbors[l_cell][3];
    }
#endif

    m_neighborKernel.computeNeighborsIntegral( m_cellInformation[l_cell].faceTypes,
                                               m_cellInformation[l_cell].faceRelations,
                                               m_cells->drMapping[l_cell],
                                               l_globalData,
                                               &m_cellData->neighboringIntegration[l_cell],
                                               l_timeIntegrated,
#ifdef ENABLE_MATRIX_PREFETCH
                                               l_faceNeighbors_prefetch,
#endif
                                               m_cells->dofs[l_cell]);
  }

#ifdef _OPENMP
  }
#endif
}

void computeDynRupGodunovState()
{
  seissol::initializers::Layer& layerData = m_dynRupTree.child(0).child<Interior>();
  DRFaceInformation*                    faceInformation                                                   = layerData.var(m_dynRup.faceInformation);
  DRGodunovData*                        godunovData                                                       = layerData.var(m_dynRup.godunovData);
  real**                                timeDerivativePlus                                                = layerData.var(m_dynRup.timeDerivativePlus);
  real**                                timeDerivativeMinus                                               = layerData.var(m_dynRup.timeDerivativeMinus);
  real                                (*godunov)[CONVERGENCE_ORDER][seissol::model::godunovState::reals]  = layerData.var(m_dynRup.godunov);

#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned face = 0; face < layerData.getNumberOfCells(); ++face) {
    unsigned prefetchFace = (face < layerData.getNumberOfCells()-1) ? face+1 : face;
    m_dynRupKernel.computeGodunovState( faceInformation[face],
                                        m_globalData,
                                       &godunovData[face],
                                        timeDerivativePlus[face],
                                        timeDerivativeMinus[face],
                                        godunov[face],
                                        timeDerivativePlus[prefetchFace],
                                        timeDerivativeMinus[prefetchFace] );
  }
}

