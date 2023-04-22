#include <Parallel/MPI.h>
#include <Solver/time_stepping/GenericGhostTimeCluster.h>

namespace seissol::time_stepping {
void GenericGhostTimeCluster::sendCopyLayer() {
  SCOREP_USER_REGION( "sendCopyLayer", SCOREP_USER_REGION_TYPE_FUNCTION )
  assert(ct.correctionTime > lastSendTime);
  lastSendTime = ct.correctionTime;
  for (unsigned int region = 0; region < meshStructure->numberOfRegions; ++region) {
    if (meshStructure->neighboringClusters[region][1] == static_cast<int>(otherGlobalClusterId)) {
     MPI_Isend(meshStructure->copyRegions[region],
                static_cast<int>(meshStructure->copyRegionSizes[region]),
                MPI_C_REAL,
                meshStructure->neighboringClusters[region][0],
                timeData + meshStructure->sendIdentifiers[region],
                seissol::MPI::mpi.comm(),
                meshStructure->sendRequests + region
               );
      sendQueue.push_back(region);
    }
  }
}

void GenericGhostTimeCluster::receiveGhostLayer() {
  SCOREP_USER_REGION( "receiveGhostLayer", SCOREP_USER_REGION_TYPE_FUNCTION )
  assert(ct.predictionTime > lastSendTime);
  for (unsigned int region = 0; region < meshStructure->numberOfRegions; ++region) {
    if (meshStructure->neighboringClusters[region][1] == static_cast<int>(otherGlobalClusterId) ) {
      MPI_Irecv(meshStructure->ghostRegions[region],
                static_cast<int>(meshStructure->ghostRegionSizes[region]),
                MPI_C_REAL,
                meshStructure->neighboringClusters[region][0],
                timeData + meshStructure->receiveIdentifiers[region],
                seissol::MPI::mpi.comm(),
                meshStructure->receiveRequests + region);
      receiveQueue.push_back(region);
    }
  }
}

bool GenericGhostTimeCluster::testQueue(MPI_Request* requests,
                                        std::list<unsigned int>& regions) {
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
}

bool GenericGhostTimeCluster::testForGhostLayerReceives() {
  SCOREP_USER_REGION( "testForGhostLayerReceives", SCOREP_USER_REGION_TYPE_FUNCTION )
  return testQueue(meshStructure->receiveRequests, receiveQueue);
}

bool GenericGhostTimeCluster::testForCopyLayerSends() {
  SCOREP_USER_REGION( "testForCopyLayerSends", SCOREP_USER_REGION_TYPE_FUNCTION )
  return testQueue(meshStructure->sendRequests, sendQueue);
}

ActResult GenericGhostTimeCluster::act() {
  // Always check for receives/send for quicker MPI progression.
  testForGhostLayerReceives();
  testForCopyLayerSends();
  return AbstractTimeCluster::act();
}

void GenericGhostTimeCluster::start() {
  assert(testForGhostLayerReceives());
  receiveGhostLayer();
}

void GenericGhostTimeCluster::predict() {
    // Doesn't do anything
}

void GenericGhostTimeCluster::correct() {
    // Doesn't do anything
}

bool GenericGhostTimeCluster::mayCorrect() {
  return testForCopyLayerSends() && AbstractTimeCluster::mayCorrect();
}

bool GenericGhostTimeCluster::mayPredict() {
  return testForGhostLayerReceives() && AbstractTimeCluster::mayPredict();
}

bool GenericGhostTimeCluster::maySync() {
  return testForGhostLayerReceives() && testForCopyLayerSends() && AbstractTimeCluster::maySync();
}

void GenericGhostTimeCluster::handleAdvancedPredictionTimeMessage(const NeighborCluster&) {
  assert(testForCopyLayerSends());
  sendCopyLayer();
}

void GenericGhostTimeCluster::handleAdvancedCorrectionTimeMessage(const NeighborCluster& neighborCluster) {
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

GenericGhostTimeCluster::GenericGhostTimeCluster(double maxTimeStepSize,
                                                 int timeStepRate,
                                                 int globalTimeClusterId,
                                                 int otherGlobalTimeClusterId,
                                                 const MeshStructure *meshStructure)
    : AbstractTimeCluster(maxTimeStepSize, timeStepRate),
      globalClusterId(globalTimeClusterId),
      otherGlobalClusterId(otherGlobalTimeClusterId),
      meshStructure(meshStructure) {}

void GenericGhostTimeCluster::reset() {
  AbstractTimeCluster::reset();
  assert(testForGhostLayerReceives());
  lastSendTime = -1;
}

void GenericGhostTimeCluster::printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) {
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
