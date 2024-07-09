#include "Parallel/MPI.h"
#include "Solver/Clustering/Communication/AbstractGhostTimeCluster.h"
#include <Common/Executor.hpp>
#include <Initializer/preProcessorMacros.hpp>
#include <Kernels/common.hpp>
#include <Solver/Clustering/AbstractTimeCluster.h>
#include <Solver/Clustering/ActorState.h>
#include <Solver/Clustering/Communication/CCLCluster.hpp>
#include <cassert>
#include <chrono>
#include <list>
#include <mpi.h>
#include <utils/logger.h>

#ifdef USE_CCL
#ifdef SEISSOL_KERNELS_CUDA
#include <nccl.h>
using StreamT = cudaStream_t;
#endif
#ifdef SEISSOL_KERNELS_HIP
#include <rccl/rccl.h>
using StreamT = hipStream_t;
#endif
#include <device.h>

#if REAL_SIZE == 8
constexpr ncclDataType_t CclReal = ncclDataType_t::ncclFloat64;
#elif REAL_SIZE == 4
constexpr ncclDataType_t CclReal = ncclDataType_t::ncclFloat32;
#endif
#endif

namespace seissol::time_stepping {

void CCLCluster::sendCopyLayer() {
#ifdef USE_CCL
  ncclGroupStart();
  for (std::size_t i = 0; i < copyClusters.size(); ++i) {
    const auto& cluster = copyClusters.at(i);
    ncclSend(cluster.data,
             cluster.size,
             CclReal, // TODO(David): use cluster.datatype here instead
             cluster.rank,
             static_cast<ncclComm_t>(cclSendComm),
             static_cast<StreamT>(streamRuntime.stream()));
  }
  ncclGroupEnd();
#endif
}

void CCLCluster::receiveGhostLayer() {
#ifdef USE_CCL
  ncclGroupStart();
  for (std::size_t i = 0; i < ghostClusters.size(); ++i) {
    const auto& cluster = ghostClusters.at(i);
    ncclRecv(cluster.data,
             cluster.size,
             CclReal, // TODO(David): use cluster.datatype here instead
             cluster.rank,
             static_cast<ncclComm_t>(cclRecvComm),
             static_cast<StreamT>(streamRuntime.stream()));
  }
  ncclGroupEnd();
#endif
}

void CCLCluster::start() { receiveGhostLayer(); }

void CCLCluster::runCompute(ComputeStep step) {
  if (step == ComputeStep::Predict) {
    streamRuntime.wait();
    sendCopyLayer();
    assert(testForCopyLayerSends());
  }
  if (step == ComputeStep::Correct) {
    streamRuntime.wait();
    assert(testForGhostLayerReceives());

    auto upcomingCorrectionSteps = ct.computeSinceLastSync.at(ComputeStep::Correct);
    if (state.step == ComputeStep::Predict && state.type == StateType::ComputeDone) {
      upcomingCorrectionSteps = ct.nextSteps();
    }

    const bool ignoreMessage = upcomingCorrectionSteps >= ct.stepsUntilSync;

    // If we are already at a sync point, we must not post an additional receive, as otherwise
    // start() posts an additional request. This is also true for the last sync point (i.e. end of
    // simulation), as in this case we do not want to have any hanging request.
    if (!ignoreMessage) {
      receiveGhostLayer();
    }
  }
}

void CCLCluster::handleAdvancedComputeTimeMessage(ComputeStep step, const NeighborCluster&) {}

CCLCluster::CCLCluster(double maxTimeStepSize,
                       int timeStepRate,
                       int globalTimeClusterId,
                       int otherGlobalTimeClusterId,
                       const MeshStructure* meshStructure,
                       void* sendComm,
                       void* recvComm)
    : CellCluster(maxTimeStepSize, timeStepRate, isDeviceOn() ? Executor::Device : Executor::Host),
      globalClusterId(globalTimeClusterId), otherGlobalClusterId(otherGlobalTimeClusterId),
      cclSendComm(sendComm), cclRecvComm(recvComm) {
  for (unsigned int region = 0; region < meshStructure->numberOfRegions; ++region) {
    if (meshStructure->neighboringClusters[region][1] == static_cast<int>(otherGlobalClusterId)) {
      copyClusters.emplace_back(RemoteCluster{
          meshStructure->copyRegions[region],
          meshStructure->copyRegionSizes[region],
          MPI_C_REAL,
          meshStructure->neighboringClusters[region][0],
          timeData + meshStructure->sendIdentifiers[region],
      });
      ghostClusters.emplace_back(RemoteCluster{
          meshStructure->ghostRegions[region],
          meshStructure->ghostRegionSizes[region],
          MPI_C_REAL,
          meshStructure->neighboringClusters[region][0],
          timeData + meshStructure->receiveIdentifiers[region],
      });

#ifdef USE_CCL
      void* sendHandle;
      ncclCommRegister(reinterpret_cast<ncclComm_t>(cclSendComm),
                       meshStructure->copyRegions[region],
                       meshStructure->copyRegionSizes[region] * sizeof(real),
                       &sendHandle);
      memorySendHandles.push_back(sendHandle);
      void* recvHandle;
      ncclCommRegister(reinterpret_cast<ncclComm_t>(cclRecvComm),
                       meshStructure->ghostRegions[region],
                       meshStructure->ghostRegionSizes[region] * sizeof(real),
                       &recvHandle);
      memoryRecvHandles.push_back(recvHandle);
#endif
    }
  }
}

void CCLCluster::reset() {
  AbstractTimeCluster::reset();
  lastSendTime = -1;
}

void CCLCluster::finalize() {
#ifdef USE_CCL
  for (void* handle : memorySendHandles) {
    ncclCommDeregister(reinterpret_cast<ncclComm_t>(cclSendComm), &handle);
  }
  for (void* handle : memoryRecvHandles) {
    ncclCommDeregister(reinterpret_cast<ncclComm_t>(cclRecvComm), &handle);
  }
#endif
}

void CCLCluster::printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) {
  const auto rank = MPI::mpi.rank();
  logError() << "Ghost: No update since " << timeSinceLastUpdate.count()
             << "[s] for global cluster " << globalClusterId << " with other cluster id "
             << otherGlobalClusterId << " at state " << actorStateToString(state)
             << " maySync = " << maySynchronize();
  for (auto& neighbor : neighbors) {
    logError() << "Neighbor with rate = " << neighbor.ct.timeStepRate
               << "PredTime = " << neighbor.ct.time.at(ComputeStep::Predict)
               << "CorrTime = " << neighbor.ct.time.at(ComputeStep::Correct)
               << "predictionsSinceSync = "
               << neighbor.ct.computeSinceLastSync.at(ComputeStep::Predict)
               << "correctionsSinceSync = "
               << neighbor.ct.computeSinceLastSync.at(ComputeStep::Correct);
  }
}

} // namespace seissol::time_stepping
