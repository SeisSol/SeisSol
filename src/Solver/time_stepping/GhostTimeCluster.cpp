#include <Parallel/MPI.h>
#include "GhostTimeCluster.h"

namespace seissol::time_stepping {
void GhostTimeCluster::sendCopyLayer(){
  SCOREP_USER_REGION( "sendCopyLayer", SCOREP_USER_REGION_TYPE_FUNCTION )
  std::cout << globalClusterId << ": sendCopyLayer, predTime = " << ct.predictionTime
  << ", corTime = " << ct.correctionTime << ", lastSendTime = " << lastSendTime
  << std::endl;
  if (ct.correctionTime <= lastSendTime) {
    std::cout << "Duplicate send!" << std::endl;
  };
  assert(ct.predictionTime > lastSendTime);
  lastSendTime = ct.correctionTime;

  for (unsigned int region = 0; region < meshStructure->numberOfRegions; ++region) {
    if (meshStructure->neighboringClusters[region][1] == static_cast<int>(globalClusterId)) {
      std::cout
      << "rank =" << seissol::MPI::mpi.rank()
      << " cluster =" << globalClusterId
      << " region =" << region << " out of " << meshStructure->numberOfRegions
      << " sending copy layer for tag " <<  timeData + meshStructure->sendIdentifiers[region]
      <<" to rank=" << meshStructure->neighboringClusters[region][0]
      << std::endl;

     MPI_Isend(meshStructure->copyRegions[region],
                static_cast<int>(meshStructure->copyRegionSizes[region]),
                MPI_C_REAL,
                meshStructure->neighboringClusters[region][0],
                timeData + meshStructure->sendIdentifiers[region],
                seissol::MPI::mpi.comm(),
                meshStructure->sendRequests + region
               );

      sendQueue.push_back(meshStructure->sendRequests + region);
    }
  }
} void GhostTimeCluster::receiveGhostLayer(){
  SCOREP_USER_REGION( "receiveGhostLayer", SCOREP_USER_REGION_TYPE_FUNCTION )
  std::cout << globalClusterId << ": ReceiveGhostLayer, predTime = " << ct.predictionTime
            << ", corTime = " << ct.correctionTime << ", lastReceiveTime = " << lastReceiveTime
            << std::endl;
  lastReceiveTime = ct.predictionTime;
  for (unsigned int region = 0; region < meshStructure->numberOfRegions; ++region) {
    if (meshStructure->neighboringClusters[region][1] == static_cast<int>(globalClusterId) ) {
      std::cout
          << "rank =" << seissol::MPI::mpi.rank()
          << " cluster =" << globalClusterId
          << " region =" << region << " out of " << meshStructure->numberOfRegions
          << " recv copy layer for tag " <<  timeData + meshStructure->receiveIdentifiers[region]
          <<" from rank=" << meshStructure->neighboringClusters[region][0]
          << std::endl;
      MPI_Irecv(meshStructure->ghostRegions[region],
                static_cast<int>(meshStructure->ghostRegionSizes[region]),
                MPI_C_REAL,
                meshStructure->neighboringClusters[region][0],
                timeData + meshStructure->receiveIdentifiers[region],
                seissol::MPI::mpi.comm(),
                meshStructure->receiveRequests + region
               );

      receiveQueue.push_back(meshStructure->receiveRequests + region );
    }
  }
}

bool GhostTimeCluster::testForGhostLayerReceives(){
  SCOREP_USER_REGION( "testForGhostLayerReceives", SCOREP_USER_REGION_TYPE_FUNCTION )

  for (auto receive = receiveQueue.begin(); receive != receiveQueue.end(); ) {
    int testSuccess = 0;
    MPI_Test(*receive, &testSuccess, MPI_STATUS_IGNORE);
    if (testSuccess) {
      std::cout << seissol::MPI::mpi.rank() << ":test for receive succ" << std::endl;
      receive = receiveQueue.erase(receive);
    } else {
      ++receive;
    }
  }
  return receiveQueue.empty();
}

bool GhostTimeCluster::testForCopyLayerSends(){
  SCOREP_USER_REGION( "testForCopyLayerSends", SCOREP_USER_REGION_TYPE_FUNCTION )

  for (auto send = sendQueue.begin(); send != sendQueue.end(); ) {
    int testSuccess = 0;
    MPI_Test(*send, &testSuccess, MPI_STATUS_IGNORE);
    if (testSuccess) {
        std::cout << seissol::MPI::mpi.rank() << ":test for send succ" << std::endl;
        send = sendQueue.erase(send);
    } else {
      ++send;
    }
  }
  return sendQueue.empty();
}

bool GhostTimeCluster::act() {
  // Always check for receives/send for quicker MPI progression.
  testForGhostLayerReceives();
  testForCopyLayerSends();
  return AbstractTimeCluster::act();
}

void GhostTimeCluster::predict() {
    // Doesn't do anything
}

void GhostTimeCluster::correct() {
    // Doesn't do anything
}
bool GhostTimeCluster::mayCorrect() {
  return testForCopyLayerSends() && AbstractTimeCluster::mayCorrect();
}
bool GhostTimeCluster::mayPredict() {
  return testForGhostLayerReceives() && AbstractTimeCluster::mayPredict();
}
void GhostTimeCluster::handleAdvancedPredictionTimeMessage(const NeighborCluster&) {
  assert(testForCopyLayerSends());
  //if(((numberOfTimeSteps + 1) % timeStepRate) == 0) {
  sendCopyLayer();
  //}
  //std::cout << "MPI recv AdvancedPredicitonTimeMessage" << std::endl;
}
void GhostTimeCluster::handleAdvancedCorrectionTimeMessage(const NeighborCluster& neighborCluster) {
  assert(testForGhostLayerReceives());
  //if(numberOfTimeSteps % timeStepRate == 0) {
  if (std::abs(neighborCluster.ct.correctionTime - syncTime) < timeTolerance) {
    std::cout << "GhostTimeCluster: ignore AdvancedCorrectionTime Message at t = "
    << neighborCluster.ct.correctionTime << ", next sync = " << syncTime << std::endl;
    return;
  } else {
  std::cout << "GhostTimeCluster: handled AdvancedCorrectionTime Message at t = "
      << neighborCluster.ct.correctionTime << ", next sync = " << syncTime << std::endl;
  }
  receiveGhostLayer();
  //}
  //std::cout << "AdvancedCorr " << ct.maxTimeStepSize << " " << msg.time << std::endl;
}

void GhostTimeCluster::cancelPendingMessages() {
  auto cancelRequest = [](MPI_Request*& request) {
    MPI_Cancel(request);
    // Note: It is neccessary to wait for cancellation and MPI_Cancel is non-blocking.
    // see: https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node72.htm
    MPI_Wait(request, MPI_STATUS_IGNORE);
  };
  assert(sendQueue.empty());
  assert(receiveQueue.empty());
  std::cout << "Begin cancelling messages" << std::endl;
  // TODO(Lukas) Should we cancel sends as well? Can be quite expensive.
  // See: https://www.mpich.org/static/docs/latest/www3/MPI_Cancel.html
  std::cout << "Cancelling sends, n=" << sendQueue.size() << std::endl;
  std::for_each(sendQueue.begin(), sendQueue.end(), cancelRequest);
  sendQueue.clear();
  std::cout << "Cancelling receives, n=" << receiveQueue.size() << std::endl;
  std::for_each(receiveQueue.begin(), receiveQueue.end(), cancelRequest);
  receiveQueue.clear();
  std::cout << "Cancelled all messages" << std::endl;


}
bool GhostTimeCluster::hasPendingMessages() {
  // & instead of && on purpose
  return !testForGhostLayerReceives() & !testForCopyLayerSends();
}
GhostTimeCluster::GhostTimeCluster(double maxTimeStepSize,
                                   int timeStepRate,
                                   double timeTolerance,
                                   int globalTimeClusterId,
                                   const MeshStructure *meshStructure)
    : AbstractTimeCluster(maxTimeStepSize, timeTolerance, timeStepRate),
      globalClusterId(globalTimeClusterId),
      meshStructure(meshStructure) {
  //reset();
  //receiveGhostLayer();
}
void GhostTimeCluster::reset() {
  std::cout << "Begin reset ghost " << globalClusterId << "\n"
  << "---------------------------" << std::endl;
  AbstractTimeCluster::reset();
  assert(testForGhostLayerReceives());
  //if (testForGhostLayerReceives()) {
    receiveGhostLayer();
  //}
  lastSendTime = -1;
  lastReceiveTime = -1;
  std::cout << "End reset ghost\n"
            << "---------------------------" << std::endl;

}

}
