#pragma once

#include "WavePropagation/dispatcher.hpp"
#include "Initializer/tree/Layer.hpp"
#include "Kernels/Time.h"
#include "Kernels/Local.h"
#include "Kernels/Neighbor.h"
#include "Kernels/Plasticity.h"
#include "Model/plasticity.hpp"
#include "implbase.hpp"

#include <vector>

namespace seissol::waveprop {
template <bool Plasticity>
class WavePropDispatcherCPU : public WavePropDispatcherPre {
  public:
  WavePropDispatcherCPU(const seissol::initializers::LTS& lts, seissol::initializers::Layer& layer)
      : WavePropDispatcherPre(lts, layer) {}

  virtual void
      dispatchPredict(double timeStepSize, double correctionTime, bool resetBuffers) override {
    // local integration buffer
    real l_integrationBuffer[tensor::I::size()] alignas(ALIGNMENT);

    // pointer for the call of the ADER-function
    real* l_bufferPointer;

    real** buffers = layer.var(lts.buffers);
    real** derivatives = layer.var(lts.derivatives);
    CellMaterialData* materialData = layer.var(lts.material);

    kernels::LocalData::Loader loader;
    loader.load(lts, layer);
    kernels::LocalTmp tmp{};

#ifdef _OPENMP
#pragma omp parallel for private(l_bufferPointer, l_integrationBuffer, tmp) schedule(static)
#endif
    for (unsigned int l_cell = 0; l_cell < layer.getNumberOfCells(); l_cell++) {
      auto data = loader.entry(l_cell);

      // We need to check, whether we can overwrite the buffer or if it is
      // needed by some other time cluster.
      // If we cannot overwrite the buffer, we compute everything in a temporary
      // local buffer and accumulate the results later in the shared buffer.
      const bool buffersProvided =
          (data.cellInformation.ltsSetup >> 8) % 2 == 1; // buffers are provided
      const bool resetMyBuffers =
          buffersProvided &&
          ((data.cellInformation.ltsSetup >> 10) % 2 == 0 || resetBuffers); // they should be reset

      if (resetMyBuffers) {
        // assert presence of the buffer
        assert(buffers[l_cell] != nullptr);

        l_bufferPointer = buffers[l_cell];
      } else {
        // work on local buffer
        l_bufferPointer = l_integrationBuffer;
      }

      timeKernel.computeAder(timeStepSize, data, tmp, l_bufferPointer, derivatives[l_cell], true);

      // Compute local integrals (including some boundary conditions)
      CellBoundaryMapping(*boundaryMapping)[4] = layer.var(lts.boundaryMapping);
      localKernel.computeIntegral(l_bufferPointer,
                                  data,
                                  tmp,
                                  &materialData[l_cell],
                                  &boundaryMapping[l_cell],
                                  correctionTime,
                                  timeStepSize);

      for (unsigned face = 0; face < 4; ++face) {
        auto& curFaceDisplacements = data.faceDisplacements[face];
        // Note: Displacement for freeSurfaceGravity is computed in Time.cpp
        if (curFaceDisplacements != nullptr &&
            data.cellInformation.faceTypes[face] != FaceType::freeSurfaceGravity) {
          kernel::addVelocity addVelocityKrnl;

          addVelocityKrnl.V3mTo2nFace = globalData.onHost->V3mTo2nFace;
          addVelocityKrnl.selectVelocity = init::selectVelocity::Values;
          addVelocityKrnl.faceDisplacement = data.faceDisplacements[face];
          addVelocityKrnl.I = l_bufferPointer;
          addVelocityKrnl.execute(face);
        }
      }

      // TODO: Integrate this step into the kernel
      // We've used a temporary buffer -> need to accumulate update in
      // shared buffer.
      if (!resetMyBuffers && buffersProvided) {
        assert(buffers[l_cell] != nullptr);

        for (unsigned int l_dof = 0; l_dof < tensor::I::size(); ++l_dof) {
          buffers[l_cell][l_dof] += l_integrationBuffer[l_dof];
        }
      }
    }
  }

  virtual void dispatchNeighborCorrect(double timeStepSize, double subTimeStart) override {
    if (layer.getNumberOfCells() == 0)
      return;

    real*(*faceNeighbors)[4] = layer.var(lts.faceNeighbors);
    CellDRMapping(*drMapping)[4] = layer.var(lts.drMapping);
    CellLocalInformation* cellInformation = layer.var(lts.cellInformation);
    auto* plasticity = layer.var(lts.plasticity);
    auto* pstrain = layer.var(lts.pstrain);
    unsigned numberOTetsWithPlasticYielding = 0;

    kernels::NeighborData::Loader loader;
    loader.load(lts, layer);

    real* l_timeIntegrated[4];
    real* l_faceNeighbors_prefetch[4];

#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(l_timeIntegrated, l_faceNeighbors_prefetch)
#endif
    for (unsigned int l_cell = 0; l_cell < layer.getNumberOfCells(); l_cell++) {
      auto data = loader.entry(l_cell);
      seissol::kernels::TimeCommon::computeIntegrals(
          timeKernel,
          data.cellInformation.ltsSetup,
          data.cellInformation.faceTypes,
          subTimeStart,
          timeStepSize,
          faceNeighbors[l_cell],
#ifdef _OPENMP
          *reinterpret_cast<real(*)[4][tensor::I::size()]>(
              &(globalData.onHost
                    ->integrationBufferLTS[omp_get_thread_num() * 4 * tensor::I::size()])),
#else
          *reinterpret_cast<real(*)[4][tensor::I::size()]>(globalData.onHost->integrationBufferLTS),
#endif
          l_timeIntegrated);

      l_faceNeighbors_prefetch[0] =
          (cellInformation[l_cell].faceTypes[1] != FaceType::dynamicRupture)
              ? faceNeighbors[l_cell][1]
              : drMapping[l_cell][1].godunov;
      l_faceNeighbors_prefetch[1] =
          (cellInformation[l_cell].faceTypes[2] != FaceType::dynamicRupture)
              ? faceNeighbors[l_cell][2]
              : drMapping[l_cell][2].godunov;
      l_faceNeighbors_prefetch[2] =
          (cellInformation[l_cell].faceTypes[3] != FaceType::dynamicRupture)
              ? faceNeighbors[l_cell][3]
              : drMapping[l_cell][3].godunov;

      // fourth face's prefetches
      if (l_cell < (layer.getNumberOfCells() - 1)) {
        l_faceNeighbors_prefetch[3] =
            (cellInformation[l_cell + 1].faceTypes[0] != FaceType::dynamicRupture)
                ? faceNeighbors[l_cell + 1][0]
                : drMapping[l_cell + 1][0].godunov;
      } else {
        l_faceNeighbors_prefetch[3] = faceNeighbors[l_cell][3];
      }

      neighborKernel.computeNeighborsIntegral(
          data, drMapping[l_cell], l_timeIntegrated, l_faceNeighbors_prefetch);

      if constexpr (Plasticity) {
        auto oneMinusIntegratingFactor = (tv > 0.0) ? 1.0 - exp(-timeStepSize / tv) : 1.0;
        numberOTetsWithPlasticYielding +=
            seissol::kernels::Plasticity::computePlasticity(oneMinusIntegratingFactor,
                                                            timeStepSize,
                                                            tv,
                                                            globalData.onHost,
                                                            &plasticity[l_cell],
                                                            data.dofs,
                                                            pstrain[l_cell]);
      }
#ifdef INTEGRATE_QUANTITIES
      seissol::SeisSol::main.postProcessor().integrateQuantities(
          timeStepSize, layer, l_cell, dofs[l_cell]);
#endif // INTEGRATE_QUANTITIES
    }
  }
};
} // namespace seissol::waveprop
