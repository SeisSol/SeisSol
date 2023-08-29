#pragma once

#include "WavePropagation/dispatcher.hpp"
#include "Initializer/tree/Layer.hpp"
#include "Equations/Kernels.hpp"
#include "Kernels/Plasticity.h"
#include "Kernels/TimeCommon.h"
#include "Model/plasticity.hpp"
#include "Common/configtensor.hpp"
#include "implbase.hpp"

#include <vector>

namespace seissol::waveprop {

template <typename Config>
class WavePropDispatcherGPU : public WavePropDispatcherPre<Config> {
  public:
  using RealT = typename Config::RealT;
  WavePropDispatcherGPU(const seissol::initializers::LTS<Config>& lts,
                        seissol::initializers::Layer& layer)
      : WavePropDispatcherPre<Config>(lts, layer) {}

#ifdef ACL_DEVICE
  void dispatchPredict(double timeStepSize, double correctionTime, bool resetBuffers) override {
    device.api->putProfilingMark("computeLocalIntegration", device::ProfilingColors::Yellow);

    RealT*(*faceNeighbors)[4] = layer.var(lts.faceNeighbors);
    auto& dataTable = layer.getConditionalTable<inner_keys::Wp>();
    auto& materialTable = layer.getConditionalTable<inner_keys::Material>();
    auto& indicesTable = layer.getConditionalTable<inner_keys::Indices>();

    kernels::LocalData<Config>::Loader loader;
    loader.load(lts, layer);
    kernels::LocalTmp<Config> tmp;

    ComputeGraphType graphType{ComputeGraphType::LocalIntegral};
    auto computeGraphKey = initializers::GraphKey(graphType, timeStepSize, true);
    auto computeGraphHandle = layer.getDeviceComputeGraphHandle(computeGraphKey);

    if (!computeGraphHandle) {
      device.api->streamBeginCapture();

      m_timeKernel.computeBatchedAder(timeStepSize, tmp, dataTable, materialTable, true);
      assert(device.api->isCircularStreamsJoinedWithDefault() &&
             "circular streams must be joined with the default stream");

      m_localKernel.computeBatchedIntegral(
          dataTable, materialTable, indicesTable, loader, tmp, timeStepSize);
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

    m_localKernel.evaluateBatchedTimeDependentBc(
        dataTable, indicesTable, loader, correctionTime, timeStepSize);

    graphType = resetBuffers ? ComputeGraphType::AccumulatedVelocities
                             : ComputeGraphType::StreamedVelocities;
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

          typename Yateto<Config>::Kernel::gpu_addVelocity displacementKrnl;
          displacementKrnl.faceDisplacement =
              entry.get(inner_keys::Wp::Id::FaceDisplacement)->getDeviceDataPtr();
          displacementKrnl.integratedVelocities = const_cast<RealT const**>(
              entry.get(inner_keys::Wp::Id::Ivelocities)->getDeviceDataPtr());
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
              Yateto<Config>::Tensor::I::Size,
              (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
              defaultStream);
        } else {
          device.algorithms.accumulateBatchedData(
              (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr(),
              (entry.get(inner_keys::Wp::Id::Buffers))->getDeviceDataPtr(),
              Yateto<Config>::Tensor::I::Size,
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

  void dispatchNeighborCorrect(double timeStepSize, double subTimeStart) override {
    device.api->putProfilingMark("computeNeighboring", device::ProfilingColors::Red);

    const double timeStepSize = timeStepSize();
    auto& table = layer.getConditionalTable<inner_keys::Wp>();

    seissol::kernels::TimeCommon::computeBatchedIntegrals(
        timeKernel, subTimeStart, timeStepSize, table);
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

    if constexpr (Config::Plasticity) {
      seissol::model::PlasticityData<>* plasticity = layer.var(lts.plasticity);
      auto oneMinusIntegratingFactor = (tv > 0.0) ? 1.0 - exp(-timeStepSize / tv) : 1.0;
      unsigned numAdjustedDofs = seissol::kernels::Plasticity::computePlasticityBatched(
          oneMinusIntegratingFactor, timeStepSize, tv, globalData.onDevice, table, plasticity);
    }

    device.api->syncDefaultStreamWithHost();
    device.api->popLastProfilingMark();
  }
#endif
};
} // namespace seissol::waveprop
