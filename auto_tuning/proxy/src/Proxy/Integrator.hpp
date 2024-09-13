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

#include "generated_code/tensor.h"

namespace tensor = seissol::tensor;
namespace kernels = seissol::kernels;

namespace proxy::cpu {
void computeAderIntegration() {
  auto& layer = ltsTree->child(0).child<Interior>();
  unsigned nrOfCells = layer.getNumberOfCells();
  real** buffers = layer.var(lts.buffers);
  real** derivatives = layer.var(lts.derivatives);

  kernels::LocalData::Loader loader;
  loader.load(lts, layer);

#ifdef _OPENMP
#pragma omp parallel
  {
    LIKWID_MARKER_START("ader");
#endif
    kernels::LocalTmp tmp(9.81);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (unsigned int cell = 0; cell < nrOfCells; cell++) {
      auto data = loader.entry(cell);
      timeKernel.computeAder(
          seissol::miniSeisSolTimeStep, data, tmp, buffers[cell], derivatives[cell]);
    }
#ifdef _OPENMP
    LIKWID_MARKER_STOP("ader");
  }
#endif
}

void computeLocalWithoutAderIntegration() {
  auto& layer = ltsTree->child(0).child<Interior>();
  unsigned nrOfCells = layer.getNumberOfCells();
  real** buffers = layer.var(lts.buffers);

  kernels::LocalData::Loader loader;
  loader.load(lts, layer);

#ifdef _OPENMP
#pragma omp parallel
  {
    LIKWID_MARKER_START("localwoader");
#endif
    kernels::LocalTmp tmp(9.81);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (unsigned int cell = 0; cell < nrOfCells; cell++) {
      auto data = loader.entry(cell);
      localKernel.computeIntegral(buffers[cell], data, tmp, nullptr, nullptr, 0, 0);
    }
#ifdef _OPENMP
    LIKWID_MARKER_STOP("localwoader");
  }
#endif
}

void computeLocalIntegration() {
  auto& layer = ltsTree->child(0).child<Interior>();
  unsigned nrOfCells = layer.getNumberOfCells();
  real** buffers = layer.var(lts.buffers);
  real** derivatives = layer.var(lts.derivatives);

  kernels::LocalData::Loader loader;
  loader.load(lts, layer);

#ifdef _OPENMP
#pragma omp parallel
  {
    LIKWID_MARKER_START("local");
#endif
    kernels::LocalTmp tmp(9.81);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (unsigned int cell = 0; cell < nrOfCells; cell++) {
      auto data = loader.entry(cell);
      timeKernel.computeAder(
          (double)seissol::miniSeisSolTimeStep, data, tmp, buffers[cell], derivatives[cell]);
      localKernel.computeIntegral(buffers[cell], data, tmp, nullptr, nullptr, 0, 0);
    }
#ifdef _OPENMP
    LIKWID_MARKER_STOP("local");
  }
#endif
}

void computeNeighboringIntegration() {
  auto& layer = ltsTree->child(0).child<Interior>();
  unsigned nrOfCells = layer.getNumberOfCells();
  real*(*faceNeighbors)[4] = layer.var(lts.faceNeighbors);
  CellDRMapping(*drMapping)[4] = layer.var(lts.drMapping);
  CellLocalInformation* cellInformation = layer.var(lts.cellInformation);

  kernels::NeighborData::Loader loader;
  loader.load(lts, layer);

  real* timeIntegrated[4];
  real* faceNeighbors_prefetch[4];

#ifdef _OPENMP
#pragma omp parallel private(timeIntegrated, faceNeighbors_prefetch)
  {
    LIKWID_MARKER_START("neighboring");
#pragma omp for schedule(static)
#endif
    for (unsigned cell = 0; cell < nrOfCells; cell++) {
      auto data = loader.entry(cell);
      seissol::kernels::TimeCommon::computeIntegrals(
          timeKernel,
          cellInformation[cell].ltsSetup,
          cellInformation[cell].faceTypes,
          0.0,
          (double)seissol::miniSeisSolTimeStep,
          faceNeighbors[cell],
#ifdef _OPENMP
          *reinterpret_cast<real(*)[4][tensor::I::size()]>(&(
              globalDataOnHost.integrationBufferLTS[omp_get_thread_num() * 4 * tensor::I::size()])),
#else
        *reinterpret_cast<real(*)[4][tensor::I::size()]>(globalDataOnHost.integrationBufferLTS),
#endif
          timeIntegrated);

      faceNeighbors_prefetch[0] = (cellInformation[cell].faceTypes[1] != FaceType::DynamicRupture)
                                      ? faceNeighbors[cell][1]
                                      : drMapping[cell][1].godunov;
      faceNeighbors_prefetch[1] = (cellInformation[cell].faceTypes[2] != FaceType::DynamicRupture)
                                      ? faceNeighbors[cell][2]
                                      : drMapping[cell][2].godunov;
      faceNeighbors_prefetch[2] = (cellInformation[cell].faceTypes[3] != FaceType::DynamicRupture)
                                      ? faceNeighbors[cell][3]
                                      : drMapping[cell][3].godunov;

      // fourth face's prefetches
      if (cell < (nrOfCells - 1)) {
        faceNeighbors_prefetch[3] =
            (cellInformation[cell + 1].faceTypes[0] != FaceType::DynamicRupture)
                ? faceNeighbors[cell + 1][0]
                : drMapping[cell + 1][0].godunov;
      } else {
        faceNeighbors_prefetch[3] = faceNeighbors[cell][3];
      }

      neighborKernel.computeNeighborsIntegral(
          data, drMapping[cell], timeIntegrated, faceNeighbors_prefetch);
    }

#ifdef _OPENMP
    LIKWID_MARKER_STOP("neighboring");
  }
#endif
}

void computeDynRupGodunovState() {
  seissol::initializer::Layer& layerData = dynRupTree->child(0).child<Interior>();
  DRFaceInformation* faceInformation = layerData.var(dynRup.faceInformation);
  DRGodunovData* godunovData = layerData.var(dynRup.godunovData);
  DREnergyOutput* drEnergyOutput = layerData.var(dynRup.drEnergyOutput);
  real** timeDerivativePlus = layerData.var(dynRup.timeDerivativePlus);
  real** timeDerivativeMinus = layerData.var(dynRup.timeDerivativeMinus);
  alignas(Alignment) real QInterpolatedPlus[ConvergenceOrder][tensor::QInterpolated::size()];
  alignas(Alignment) real QInterpolatedMinus[ConvergenceOrder][tensor::QInterpolated::size()];

#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(QInterpolatedPlus, QInterpolatedMinus)
#endif
  for (unsigned face = 0; face < layerData.getNumberOfCells(); ++face) {
    unsigned prefetchFace = (face < layerData.getNumberOfCells() - 1) ? face + 1 : face;
    dynRupKernel.spaceTimeInterpolation(faceInformation[face],
                                        &globalDataOnHost,
                                        &godunovData[face],
                                        &drEnergyOutput[face],
                                        timeDerivativePlus[face],
                                        timeDerivativeMinus[face],
                                        QInterpolatedPlus,
                                        QInterpolatedMinus,
                                        timeDerivativePlus[prefetchFace],
                                        timeDerivativeMinus[prefetchFace]);
  }
}
} // namespace proxy::cpu
