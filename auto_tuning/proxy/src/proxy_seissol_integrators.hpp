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

void registerMarkers() {
    #pragma omp parallel
    {
        LIKWID_MARKER_REGISTER("ader");
        LIKWID_MARKER_REGISTER("localwoader");
        LIKWID_MARKER_REGISTER("local");
        LIKWID_MARKER_REGISTER("neighboring");
    }
}

namespace proxy::cpu {
  void computeAderIntegration() {
    auto&                 layer           = m_ltsTree->child(0).child<Interior>();
    unsigned              nrOfCells       = layer.getNumberOfCells();
    real**                buffers                       = layer.var(m_lts.buffers);
    real**                derivatives                   = layer.var(m_lts.derivatives);

    kernels::LocalData::Loader loader;
    loader.load(m_lts, layer);

  #ifdef _OPENMP
    #pragma omp parallel
    {
    LIKWID_MARKER_START("ader");
    kernels::LocalTmp tmp;
    #pragma omp for schedule(static)
  #endif
    for( unsigned int l_cell = 0; l_cell < nrOfCells; l_cell++ ) {
      auto data = loader.entry(l_cell);
      m_timeKernel.computeAder(              seissol::miniSeisSolTimeStep,
                                             data,
                                             tmp,
                                             buffers[l_cell],
                                             derivatives[l_cell] );
    }
  #ifdef _OPENMP
    LIKWID_MARKER_STOP("ader");
    }
  #endif
  }

  void computeLocalWithoutAderIntegration() {
    auto&                 layer           = m_ltsTree->child(0).child<Interior>();
    unsigned              nrOfCells       = layer.getNumberOfCells();
    real**                buffers                       = layer.var(m_lts.buffers);

    kernels::LocalData::Loader loader;
    loader.load(m_lts, layer);

  #ifdef _OPENMP
    #pragma omp parallel
    {
    LIKWID_MARKER_START("localwoader");
    kernels::LocalTmp tmp;
  #pragma omp for schedule(static)
  //#pragma omp taskloop grainsize(4) untied
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
    LIKWID_MARKER_STOP("localwoader");
    }
  #endif
  }

  template<typename Function>
  void scheduleTasksRecursively(Function&& f, int from, int to, int level = 0) {
    const auto nr = to - from;
    const auto chunkSize = 5;
    if (nr > 1 && (1 << (level - 1)) <= omp_get_num_threads()) {
      const auto half = from + nr/2;
#pragma omp task firstprivate(from, half, level) untied
      scheduleTasksRecursively(f, from, half, level + 1);
      scheduleTasksRecursively(f, half, to, level + 1);
    } else {
        f(from, to);
    }
  }

  __attribute__((always_inline)) inline void computeLocalIntegrationKernel(kernels::LocalData::Loader& loader,
                                     kernels::LocalTmp& tmp,
                                     real** buffers,
                                     real** derivatives,
                                     unsigned cell) {
    auto data = loader.entry(cell);
    m_timeKernel.computeAder((double)seissol::miniSeisSolTimeStep,
                             data,
                             tmp,
                             buffers[cell],
                             derivatives[cell] );
    m_localKernel.computeIntegral(buffers[cell],
                                  data,
                                  tmp,
                                  nullptr,
                                  nullptr,
                                  0,
                                  0);
  }

  void computeLocalIntegration(const ProxyKernelConfig& config) {
    auto&                 layer           = m_ltsTree->child(0).child<Interior>();
    unsigned              nrOfCells       = layer.getNumberOfCells();
    real**                buffers                       = layer.var(m_lts.buffers);
    real**                derivatives                   = layer.var(m_lts.derivatives);

    kernels::LocalData::Loader loader;
    loader.load(m_lts, layer);

#pragma omp parallel
    {
      LIKWID_MARKER_START("local");
      kernels::LocalTmp tmp;
      switch (config.parallelizationStrategy) {
      case ParallelizationStrategy::ParallelFor:
#pragma omp for schedule(static)
      for (unsigned int cell = 0; cell < nrOfCells; cell++) {
        computeLocalIntegrationKernel(loader, tmp, buffers, derivatives, cell);
      }
        break;
      case ParallelizationStrategy::Taskloop:
#pragma omp single
#pragma omp taskloop untied grainsize(config.grainsize)
        for (unsigned int cell = 0; cell < nrOfCells; cell++) {
          computeLocalIntegrationKernel(loader, tmp, buffers, derivatives, cell);
        }
        break;
      }
      LIKWID_MARKER_STOP("local");
    }
  }

__attribute__((always_inline)) inline void computeNeighboringIntegrationKernel(kernels::NeighborData::Loader loader,
                                          CellLocalInformation* cellInformation,
                                          real* (*faceNeighbors)[4],
                                          CellDRMapping (*drMapping)[4],
                                          real* l_timeIntegrated[4],
                                          real* l_faceNeighbors_prefetch[4],
                                          unsigned nrOfCells,
                                          unsigned cell) {
    auto data = loader.entry(cell);
    seissol::kernels::TimeCommon::computeIntegrals(
        m_timeKernel,
        cellInformation[cell].ltsSetup,
        cellInformation[cell].faceTypes,
        0.0,
        (double)seissol::miniSeisSolTimeStep,
        faceNeighbors[cell],
        *reinterpret_cast<real(*)[4][tensor::I::size()]>(&(
            m_globalDataOnHost.integrationBufferLTS[omp_get_thread_num() * 4 * tensor::I::size()])),
        l_timeIntegrated);

    l_faceNeighbors_prefetch[0] = (cellInformation[cell].faceTypes[1] != FaceType::dynamicRupture)
                                      ? faceNeighbors[cell][1]
                                      : drMapping[cell][1].godunov;
    l_faceNeighbors_prefetch[1] = (cellInformation[cell].faceTypes[2] != FaceType::dynamicRupture)
                                      ? faceNeighbors[cell][2]
                                      : drMapping[cell][2].godunov;
    l_faceNeighbors_prefetch[2] = (cellInformation[cell].faceTypes[3] != FaceType::dynamicRupture)
                                      ? faceNeighbors[cell][3]
                                      : drMapping[cell][3].godunov;

    // fourth face's prefetches
    if (cell < (nrOfCells - 1)) {
      l_faceNeighbors_prefetch[3] =
          (cellInformation[cell + 1].faceTypes[0] != FaceType::dynamicRupture)
              ? faceNeighbors[cell + 1][0]
              : drMapping[cell + 1][0].godunov;
    } else {
      l_faceNeighbors_prefetch[3] = faceNeighbors[cell][3];
    }

    m_neighborKernel.computeNeighborsIntegral(
        data, drMapping[cell], l_timeIntegrated, l_faceNeighbors_prefetch);
  }

  void computeNeighboringIntegration(const ProxyKernelConfig& config) {
    auto& layer = m_ltsTree->child(0).child<Interior>();
    unsigned nrOfCells = layer.getNumberOfCells();
    real*(*faceNeighbors)[4] = layer.var(m_lts.faceNeighbors);
    CellDRMapping(*drMapping)[4] = layer.var(m_lts.drMapping);
    CellLocalInformation* cellInformation = layer.var(m_lts.cellInformation);

    kernels::NeighborData::Loader loader;
    loader.load(m_lts, layer);

    real* l_timeIntegrated[4];
    real* l_faceNeighbors_prefetch[4];

#pragma omp parallel private(l_timeIntegrated, l_faceNeighbors_prefetch)
    {
      LIKWID_MARKER_START("neighboring");

      switch (config.parallelizationStrategy) {
      case ParallelizationStrategy::ParallelFor:
#pragma omp for schedule(static)
        for (unsigned cell = 0; cell < nrOfCells; cell++) {
          computeNeighboringIntegrationKernel(loader,
                                              cellInformation,
                                              faceNeighbors,
                                              drMapping,
                                              l_timeIntegrated,
                                              l_faceNeighbors_prefetch,
                                              nrOfCells,
                                              cell);
        }
        break;
      case ParallelizationStrategy::Taskloop:
#pragma omp single
#pragma omp taskloop untied grainsize(config.grainsize)
        for (unsigned cell = 0; cell < nrOfCells; cell++) {
          computeNeighboringIntegrationKernel(loader,
                                              cellInformation,
                                              faceNeighbors,
                                              drMapping,
                                              l_timeIntegrated,
                                              l_faceNeighbors_prefetch,
                                              nrOfCells,
                                              cell);
        }
        break;
      }

      LIKWID_MARKER_STOP("neighboring");
    }
  }


  void computeDynRupGodunovState()
  {
    seissol::initializers::Layer& layerData = m_dynRupTree->child(0).child<Interior>();
    DRFaceInformation* faceInformation = layerData.var(m_dynRup.faceInformation);
    DRGodunovData* godunovData = layerData.var(m_dynRup.godunovData);
    DREnergyOutput* drEnergyOutput = layerData.var(m_dynRup.drEnergyOutput);
    real** timeDerivativePlus = layerData.var(m_dynRup.timeDerivativePlus);
    real** timeDerivativeMinus = layerData.var(m_dynRup.timeDerivativeMinus);
    alignas(ALIGNMENT) real QInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()];
    alignas(ALIGNMENT) real QInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()];

  #ifdef _OPENMP
    #pragma omp parallel for schedule(static) private(QInterpolatedPlus,QInterpolatedMinus)
  #endif
    for (unsigned face = 0; face < layerData.getNumberOfCells(); ++face) {
      unsigned prefetchFace = (face < layerData.getNumberOfCells()-1) ? face+1 : face;
      m_dynRupKernel.spaceTimeInterpolation(  faceInformation[face],
                                             &m_globalDataOnHost,
                                             &godunovData[face],
                                             &drEnergyOutput[face],
                                              timeDerivativePlus[face],
                                              timeDerivativeMinus[face],
                                              QInterpolatedPlus,
                                              QInterpolatedMinus,
                                              timeDerivativePlus[prefetchFace],
                                              timeDerivativeMinus[prefetchFace] );
    }
  }
} // namespace proxy::cpu
