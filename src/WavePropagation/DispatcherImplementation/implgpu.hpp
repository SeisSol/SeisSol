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
    template<bool Plasticity>
    class WavePropDispatcherGPU : public WavePropDispatcherPre {
    public:
      WavePropDispatcherGPU(const seissol::initializers::LTS& lts, seissol::initializers::Layer& layer)
            : WavePropDispatcherPre(lts, layer) {}

#ifdef ACL_DEVICE
        virtual void dispatchPredict(double timeStepSize, double correctionTime, bool resetBuffers) override {
            device.api->putProfilingMark("computeLocalIntegration", device::ProfilingColors::Yellow);

  real* (*faceNeighbors)[4] = layer.var(lts.faceNeighbors);
  auto& dataTable = layer.getConditionalTable<inner_keys::Wp>();
  auto& materialTable = layer.getConditionalTable<inner_keys::Material>();
  auto& indicesTable = layer.getConditionalTable<inner_keys::Indices>();

  kernels::LocalData::Loader loader;
  loader.load(lts, layer);
  kernels::LocalTmp tmp;

  ComputeGraphType graphType{ComputeGraphType::LocalIntegral};
  auto computeGraphKey = initializers::GraphKey(graphType, timeStepSize, true);
  auto computeGraphHandle = layer.getDeviceComputeGraphHandle(computeGraphKey);

  if (!computeGraphHandle) {
    device.api->streamBeginCapture();

    m_timeKernel.computeBatchedAder(timeStepSize,
                                    tmp,
                                    dataTable,
                                    materialTable,
                                    true);
    assert(device.api->isCircularStreamsJoinedWithDefault() &&
           "circular streams must be joined with the default stream");

    m_localKernel.computeBatchedIntegral(dataTable,
                                         materialTable,
                                         indicesTable,
                                         loader,
                                         tmp,
                                         timeStepSize);
    assert(device.api->isCircularStreamsJoinedWithDefault() &&
           "circular streams must be joined with the default stream");

    device.api->streamEndCapture();

    computeGraphHandle = device.api->getLastGraphHandle();
    layer.updateDeviceComputeGraphHandle(computeGraphKey, computeGraphHandle);
    device.api->syncDefaultStreamWithHost();
  }

  if (computeGraphHandle.isInitialized()) {
    device.api->launchGraph(computeGraphHandle);
    device.api->syncGraph(computeGraphHandle);
  }

  m_localKernel.evaluateBatchedTimeDependentBc(dataTable,
                                               indicesTable,
                                               loader,
                                               correctionTime,
                                               timeStepSize);

  graphType = resetBuffers ? ComputeGraphType::AccumulatedVelocities : ComputeGraphType::StreamedVelocities;
  computeGraphKey = initializers::GraphKey(graphType);
  computeGraphHandle = layer.getDeviceComputeGraphHandle(computeGraphKey);

  if (!computeGraphHandle) {
    device.api->streamBeginCapture();

    auto defaultStream = device.api->getDefaultStream();

    for (unsigned face = 0; face < 4; ++face) {
      ConditionalKey key(*KernelNames::FaceDisplacements, *ComputationKind::None, face);
      if (dataTable.find(key) != dataTable.end()) {
        auto& entry = dataTable[key];
        // NOTE: integrated velocities have been computed implicitly, i.e
        // it is 6th, 7the and 8th columns of integrated dofs

        kernel::gpu_addVelocity displacementKrnl;
        displacementKrnl.faceDisplacement = entry.get(inner_keys::Wp::Id::FaceDisplacement)->getDeviceDataPtr();
        displacementKrnl.integratedVelocities = const_cast<real const**>(entry.get(inner_keys::Wp::Id::Ivelocities)->getDeviceDataPtr());
        displacementKrnl.V3mTo2nFace = globalData.onDevice->V3mTo2nFace;

        // Note: this kernel doesn't require tmp. memory
        displacementKrnl.numElements = entry.get(inner_keys::Wp::Id::FaceDisplacement)->getSize();
        displacementKrnl.streamPtr = defaultStream;
        displacementKrnl.execute(face);
      }
    }

    ConditionalKey key = ConditionalKey(*KernelNames::Time, *ComputationKind::WithLtsBuffers);
    if (dataTable.find(key) != dataTable.end()) {
      auto& entry = dataTable[key];

      if (resetBuffers) {
        device.algorithms.streamBatchedData(
            (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr(),
            (entry.get(inner_keys::Wp::Id::Buffers))->getDeviceDataPtr(),
            tensor::I::Size,
            (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
            defaultStream);
      } else {
        device.algorithms.accumulateBatchedData(
            (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr(),
            (entry.get(inner_keys::Wp::Id::Buffers))->getDeviceDataPtr(),
            tensor::I::Size,
            (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
            defaultStream);
      }
    }

    device.api->streamEndCapture();

    computeGraphHandle = device.api->getLastGraphHandle();
    layer.updateDeviceComputeGraphHandle(computeGraphKey, computeGraphHandle);
    device.api->syncDefaultStreamWithHost();
  }

  if (computeGraphHandle.isInitialized()) {
    device.api->launchGraph(computeGraphHandle);
    device.api->syncGraph(computeGraphHandle);
  }

  device.api->popLastProfilingMark();
        }

        virtual void dispatchNeighborCorrect(double timeStepSize, double subTimeStart) override {
          device.api->putProfilingMark("computeNeighboring", device::ProfilingColors::Red);

          const double timeStepSize = timeStepSize();
          auto& table = layer.getConditionalTable<inner_keys::Wp>();

          seissol::kernels::TimeCommon::computeBatchedIntegrals(timeKernel,
                                                                subTimeStart,
                                                                timeStepSize,
                                                                table);
          assert(device.api->isCircularStreamsJoinedWithDefault() &&
                "circular streams must be joined with the default stream");

          ComputeGraphType graphType = ComputeGraphType::NeighborIntegral;
          auto computeGraphKey = initializers::GraphKey(graphType);
          auto computeGraphHandle = layer.getDeviceComputeGraphHandle(computeGraphKey);

          if (!computeGraphHandle) {
            device.api->streamBeginCapture();

            neighborKernel.computeBatchedNeighborsIntegral(table);
            assert(device.api->isCircularStreamsJoinedWithDefault() &&
                  "circular streams must be joined with the default stream");

            device.api->streamEndCapture();

            computeGraphHandle = device.api->getLastGraphHandle();
            layer.updateDeviceComputeGraphHandle(computeGraphKey, computeGraphHandle);
            device.api->syncDefaultStreamWithHost();
          }

          if (computeGraphHandle.isInitialized()) {
            // Note: graph stream needs to wait the default stream
            // (used in `computeBatchedIntegrals`)
            device.api->syncDefaultStreamWithHost();
            device.api->launchGraph(computeGraphHandle);
            device.api->syncGraph(computeGraphHandle);
          }

          if constexpr (Plasticity) {
            seissol::model::PlasticityData<>* plasticity = layer.var(lts.plasticity);
            auto oneMinusIntegratingFactor = (tv > 0.0) ? 1.0 - exp(-timeStepSize / tv) : 1.0;
            unsigned numAdjustedDofs = seissol::kernels::Plasticity::computePlasticityBatched(oneMinusIntegratingFactor,
                                                                                              timeStepSize,
                                                                                              tv,
                                                                                              globalData.onDevice,
                                                                                              table,
                                                                                              plasticity);
          }

          device.api->syncDefaultStreamWithHost();
          device.api->popLastProfilingMark();
        }
#endif
    };
}
