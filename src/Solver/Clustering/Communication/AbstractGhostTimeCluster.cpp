#include "Solver/Clustering/Communication/AbstractGhostTimeCluster.h"
#include "Parallel/MPI.h"
#include <Common/Executor.hpp>
#include <Initializer/preProcessorMacros.hpp>
#include <Kernels/common.hpp>
#include <Solver/Clustering/AbstractTimeCluster.h>
#include <Solver/Clustering/ActorState.h>
#include <cassert>
#include <chrono>
#include <list>
#include <mpi.h>
#include <utils/logger.h>

namespace seissol::time_stepping {
bool AbstractGhostTimeCluster::testQueue(std::vector<MPI_Request>& requests,
                                         std::list<std::size_t>& remaining) {
  // we may now also use MPI_Testall here---if that's as efficient as the other method
  // (and doesn't only poll the first non-finished request)

  for (auto region = remaining.begin(); region != remaining.end();) {
    MPI_Request* request = &requests[*region];
    int testSuccess = 0;
    MPI_Test(request, &testSuccess, MPI_STATUS_IGNORE);
    if (testSuccess) {
      region = remaining.erase(region);
    } else {
      ++region;
    }
  }
  return remaining.empty();
}

bool AbstractGhostTimeCluster::testForCopyLayerSends() {
  SCOREP_USER_REGION("testForCopyLayerSends", SCOREP_USER_REGION_TYPE_FUNCTION)
  return testQueue(sendRequests, sendQueue);
}

ActResult AbstractGhostTimeCluster::act() {
  // Always check for receives/send for quicker MPI progression.
  testForGhostLayerReceives();
  testForCopyLayerSends();
  return AbstractTimeCluster::act();
}

void AbstractGhostTimeCluster::start() {
  assert(testForGhostLayerReceives());
  receiveGhostLayer();
}

void AbstractGhostTimeCluster::runCompute(ComputeStep step) {
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

bool AbstractGhostTimeCluster::pollCompute(ComputeStep step) {
  if (step == ComputeStep::Predict) {
    return testForCopyLayerSends();
  }
  if (step == ComputeStep::Correct) {
    return testForGhostLayerReceives();
  }
  return true;
}

void AbstractGhostTimeCluster::handleAdvancedComputeTimeMessage(ComputeStep step,
                                                                const NeighborCluster&) {}

AbstractGhostTimeCluster::AbstractGhostTimeCluster(double maxTimeStepSize,
                                                   int timeStepRate,
                                                   int globalTimeClusterId,
                                                   int otherGlobalTimeClusterId,
                                                   const MeshStructure* meshStructure)
    : CellCluster(maxTimeStepSize, timeStepRate, isDeviceOn() ? Executor::Device : Executor::Host),
      globalClusterId(globalTimeClusterId), otherGlobalClusterId(otherGlobalTimeClusterId),
      meshStructure(meshStructure) {
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
      sendRequests.emplace_back(MPI_REQUEST_NULL);
      recvRequests.emplace_back(MPI_REQUEST_NULL);
    }
  }
}

void AbstractGhostTimeCluster::reset() {
  AbstractTimeCluster::reset();
  assert(testForGhostLayerReceives());
  lastSendTime = -1;
}

void AbstractGhostTimeCluster::printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) {
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
