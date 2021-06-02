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

#include <generated_code/tensor.h>

namespace tensor = seissol::tensor;
namespace kernels = seissol::kernels;

namespace proxy::cpu {
  void computeAderIntegration() {
    /*
    auto&                 layer           = m_ltsTree->child(0).child<Interior>();
    unsigned              nrOfCells       = layer.getNumberOfCells();
    real**                buffers                       = layer.var(m_lts.buffers);
    real**                derivatives                   = layer.var(m_lts.derivatives);

    kernels::LocalData::Loader loader;
    loader.load(m_lts, layer);
     */

    /*
#ifdef _OPENMP
  #pragma omp parallel
  {
  kernels::LocalTmp tmp;
  #pragma omp for schedule(static)
#endif
  for( unsigned int l_cell = 0; l_cell < nrOfCells; l_cell++ ) {
    auto data = loader.entry(l_cell);
    m_timeKernel.computeAder(              m_timeStepWidthSimulation,
                                           data,
                                           tmp,
                                           buffers[l_cell],
                                           derivatives[l_cell] );
  }
     */
  #ifdef _OPENMP
    }
  #endif
  }

  void computeLocalWithoutAderIntegration() {
  /*
    auto&                 layer           = m_ltsTree->child(0).child<Interior>();
    unsigned              nrOfCells       = layer.getNumberOfCells();
    real**                buffers                       = layer.var(m_lts.buffers);

    kernels::LocalData::Loader loader;
    loader.load(m_lts, layer);
    */

    /*
#ifdef _OPENMP
  #pragma omp parallel
  {
  kernels::LocalTmp tmp;
  #pragma omp for schedule(static)
#endif
  for( unsigned int l_cell = 0; l_cell < nrOfCells; l_cell++ ) {
    auto data = loader.entry(l_cell);
    m_localKernel.computeIntegral(buffers[l_cell],
                                  data,
                                  tmp,
                                  nullptr,
                                  nullptr,
                                  0,
                                  0);
    }
  #ifdef _OPENMP
    }
  #endif
                                  */
  }

  void computeLocalIntegration() {
    auto elementViewInterior = proxyData->getElementView();
    const auto nrOfCells = elementViewInterior.size();

  #ifdef _OPENMP
    #pragma omp parallel default(none) shared(nrOfCells, elementViewInterior)
    {
    kernels::LocalTmp tmp;
    #pragma omp for schedule(static)
  #endif
    for( unsigned int l_cell = 0; l_cell < nrOfCells; l_cell++ ) {
      auto& curElement = elementViewInterior[l_cell];
      auto curDerivatives = nullptr;
      auto* curBuffers = curElement.get<buffer>();
      auto* curDofs = curElement.get<dofs>().data();
      auto* curDisplacements = curElement.get<displacements>();
      auto& curLocalIntegration = curElement.get<localIntegrationData>();
      auto& curNeighboringIntegration = curElement.get<neighborIntegrationData>();
      auto& curCellInformation = curElement.get<cellLocalInformation>();

      m_timeKernel.computeAder((double)m_timeStepWidthSimulation,
                               curDofs,
                               curDisplacements,
                               curLocalIntegration,
                               curCellInformation,
                               tmp,
                               curBuffers,
                               curDerivatives);
      m_localKernel.computeIntegral(curBuffers,
                                    curDofs,
                                    curDisplacements,
                                    curLocalIntegration,
                                    curNeighboringIntegration,
                                    curCellInformation,
                                    tmp,
                                    nullptr,
                                    nullptr,
                                    0,
                                    0);
    }
  #ifdef _OPENMP
    }
  #endif
  }

  void computeNeighboringIntegration() {
    auto elementViewInterior = proxyData->getElementView();
    const auto nrOfCells = elementViewInterior.size();

    real *l_timeIntegrated[4];
  #ifdef ENABLE_MATRIX_PREFETCH
    real *l_faceNeighbors_prefetch[4];
  #endif

  #ifdef _OPENMP
  #  ifdef ENABLE_MATRIX_PREFETCH
    #pragma omp parallel default(none) private(l_timeIntegrated, l_faceNeighbors_prefetch) shared(nrOfCells, elementViewInterior)
  #  else
    #pragma omp parallel private(l_timeIntegrated)
  #  endif
    {
    #pragma omp for schedule(static)
  #endif
    for( unsigned l_cell = 0; l_cell < nrOfCells; l_cell++ ) {
      auto& curElement = elementViewInterior[l_cell];
      auto curDerivatives = nullptr;
      auto* curBuffers = curElement.get<buffer>();
      auto* curDofs = curElement.get<dofs>().data();
      auto* curDisplacements = curElement.get<displacements>();
      auto& curLocalIntegration = curElement.get<localIntegrationData>();
      auto& curNeighboringIntegration = curElement.get<neighborIntegrationData>();
      auto& curCellInformation = curElement.get<cellLocalInformation>();
      auto& curFaceNeighbors = curElement.get<faceNeighbors>();
      auto& curDrMapping = curElement.get<cellDrMapping>();
      seissol::kernels::TimeCommon::computeIntegrals( m_timeKernel,
                                                      curCellInformation.ltsSetup,
                                                      curCellInformation.faceTypes,
                                                      0.0,
                                              (double)m_timeStepWidthSimulation,
                                                      curFaceNeighbors.data(),
  #ifdef _OPENMP
                                                      *reinterpret_cast<real (*)[4][tensor::I::size()]>(&(m_globalDataOnHost.integrationBufferLTS[omp_get_thread_num()*4*tensor::I::size()])),
  #else
                                                      *reinterpret_cast<real (*)[4][tensor::I::size()]>(m_globalData.integrationBufferLTS),
  #endif
                                                      l_timeIntegrated );

  #ifdef ENABLE_MATRIX_PREFETCH
  #pragma message("the current prefetch structure (flux matrices and tDOFs is tuned for higher order and shouldn't be harmful for lower orders")
      l_faceNeighbors_prefetch[0] = (curCellInformation.faceTypes[1] != FaceType::dynamicRupture)
          ? curFaceNeighbors[1] : curDrMapping[1].godunov;
      l_faceNeighbors_prefetch[1] = (curCellInformation.faceTypes[2] != FaceType::dynamicRupture)
          ? curFaceNeighbors[2] : curDrMapping[2].godunov;
      l_faceNeighbors_prefetch[2] = (curCellInformation.faceTypes[3] != FaceType::dynamicRupture)
          ? curFaceNeighbors[3] : curDrMapping[3].godunov;

      // fourth face's prefetches
      if (l_cell < (nrOfCells-1) ) {
        auto& nextElement = elementViewInterior[l_cell];
        auto& nextCellInformation = nextElement.get<cellLocalInformation>();
        auto& nextDrMapping = nextElement.get<cellDrMapping>();
        auto& nextFaceNeighbors = nextElement.get<faceNeighbors>();
        l_faceNeighbors_prefetch[3] = (nextCellInformation.faceTypes[0] != FaceType::dynamicRupture) ?
                                      nextFaceNeighbors[0] : nextDrMapping[0].godunov;
      } else {
        l_faceNeighbors_prefetch[3] = curFaceNeighbors[3];
      }
  #endif

      m_neighborKernel.computeNeighborsIntegral(curDofs,
                                                curNeighboringIntegration,
                                                curCellInformation,
                                                curDrMapping.data(),
                                                l_timeIntegrated,
                                                l_faceNeighbors_prefetch
                                                );
    }

  #ifdef _OPENMP
    }
  #endif

  }

  void computeDynRupGodunovState()
  {
  /*
    seissol::initializers::Layer& layerData = m_dynRupTree->child(0).child<Interior>();
    DRFaceInformation* faceInformation = layerData.var(m_dynRup.faceInformation);
    DRGodunovData* godunovData = layerData.var(m_dynRup.godunovData);
    real** timeDerivativePlus = layerData.var(m_dynRup.timeDerivativePlus);
    real** timeDerivativeMinus = layerData.var(m_dynRup.timeDerivativeMinus);
    alignas(ALIGNMENT) real QInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()];
    alignas(ALIGNMENT) real QInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()];
    */

  /*
  #ifdef _OPENMP
    #pragma omp parallel for schedule(static) private(QInterpolatedPlus,QInterpolatedMinus)
  #endif
    for (unsigned face = 0; face < layerData.getNumberOfCells(); ++face) {
      unsigned prefetchFace = (face < layerData.getNumberOfCells()-1) ? face+1 : face;
      m_dynRupKernel.spaceTimeInterpolation(  faceInformation[face],
                                             &m_globalDataOnHost,
                                             &godunovData[face],
                                              timeDerivativePlus[face],
                                              timeDerivativeMinus[face],
                                              QInterpolatedPlus,
                                              QInterpolatedMinus,
                                              timeDerivativePlus[prefetchFace],
                                              timeDerivativeMinus[prefetchFace] );
    }
  }
    */
} // namespace proxy::cpu
