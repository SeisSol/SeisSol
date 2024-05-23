#include "Parallel/MPI.h"
#include "Solver/time_stepping/AbstractGhostTimeCluster.h"

#ifdef USE_CCL
#include <device.h>
#endif


namespace seissol::time_stepping {
bool AbstractGhostTimeCluster::testQueue(MPI_Request* requests,
                                         std::list<unsigned int>& regions) {
#ifdef USE_CCL
  return true; // should never be called
#else
  for (auto region = regions.begin(); region != regions.end();) {
    MPI_Request *request = &requests[*region];
    int testSuccess = 0;
    MPI_Test(request, &testSuccess, MPI_STATUS_IGNORE);
    if (testSuccess) {
      region = regions.erase(region);
    } else {
      ++region;
    }
  }
  return regions.empty();
#endif
}

bool AbstractGhostTimeCluster::testForCopyLayerSends() {
  SCOREP_USER_REGION( "testForCopyLayerSends", SCOREP_USER_REGION_TYPE_FUNCTION )
#ifdef USE_CCL
  return ::device::DeviceInstance::getInstance().api->isStreamWorkDone(sendStream);
#else
  return testQueue(meshStructure->sendRequests, sendQueue);
#endif
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

void AbstractGhostTimeCluster::predict() {
  // Doesn't do anything
}

void AbstractGhostTimeCluster::correct() {
  // Doesn't do anything
}

bool AbstractGhostTimeCluster::mayCorrect() {
  return testForCopyLayerSends() && AbstractTimeCluster::mayCorrect();
}

bool AbstractGhostTimeCluster::mayPredict() {
  return testForGhostLayerReceives() && AbstractTimeCluster::mayPredict();
}

bool AbstractGhostTimeCluster::maySync() {
  return testForGhostLayerReceives() && testForCopyLayerSends() && AbstractTimeCluster::maySync();
}

void AbstractGhostTimeCluster::handleAdvancedPredictionTimeMessage(const NeighborCluster&) {
  assert(testForCopyLayerSends());
  sendCopyLayer();
}

void AbstractGhostTimeCluster::handleAdvancedCorrectionTimeMessage(const NeighborCluster& neighborCluster) {
  assert(testForGhostLayerReceives());

  auto upcomingCorrectionSteps = ct.stepsSinceLastSync;
  if (state == ActorState::Predicted) {
    upcomingCorrectionSteps = ct.nextCorrectionSteps();
  }

  const bool ignoreMessage = upcomingCorrectionSteps >= ct.stepsUntilSync;

  // If we are already at a sync point, we must not post an additional receive, as otherwise start() posts an additional
  // request!
  // This is also true for the last sync point (i.e. end of simulation), as in this case we do not want to have any
  // hanging request.
  if (!ignoreMessage) {
    receiveGhostLayer();
  }
}

AbstractGhostTimeCluster::AbstractGhostTimeCluster(double maxTimeStepSize,
                                                 int timeStepRate,
                                                 int globalTimeClusterId,
                                                 int otherGlobalTimeClusterId,
                                                 const MeshStructure *meshStructure)
    : AbstractTimeCluster(maxTimeStepSize, timeStepRate),
      globalClusterId(globalTimeClusterId),
      otherGlobalClusterId(otherGlobalTimeClusterId),
      meshStructure(meshStructure) {
#ifdef USE_CCL
  sendStream = ::device::DeviceInstance::getInstance().api->createGenericStream();
  recvStream = ::device::DeviceInstance::getInstance().api->createGenericStream();
#endif
}

void AbstractGhostTimeCluster::finalize() {
#ifdef USE_CCL
::device::DeviceInstance::getInstance().api->destroyGenericStream(sendStream);
::device::DeviceInstance::getInstance().api->destroyGenericStream(recvStream);
#endif
}

void AbstractGhostTimeCluster::reset() {
  AbstractTimeCluster::reset();
  assert(testForGhostLayerReceives());
  lastSendTime = -1;
}

void AbstractGhostTimeCluster::printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) {
  const auto rank = MPI::mpi.rank();
  logError()
      << "Ghost: No update since " << timeSinceLastUpdate.count()
      << "[s] for global cluster " << globalClusterId
      << " with other cluster id " << otherGlobalClusterId
      << " at state " << actorStateToString(state)
      << " mayPredict = " << mayPredict()
      << " mayPredict (steps) = " << AbstractTimeCluster::mayPredict()
      << " mayCorrect = " << mayCorrect()
      << " mayCorrect (steps) = " << AbstractTimeCluster::mayCorrect()
      << " maySync = " << maySync();
  for (auto& neighbor : neighbors) {
    logError()
        << "Neighbor with rate = " << neighbor.ct.timeStepRate
        << "PredTime = " << neighbor.ct.predictionTime
        << "CorrTime = " << neighbor.ct.correctionTime
        << "predictionsSinceSync = " << neighbor.ct.predictionsSinceLastSync
        << "correctionsSinceSync = " << neighbor.ct.stepsSinceLastSync;
  }
}

} // namespace seissol::time_stepping
