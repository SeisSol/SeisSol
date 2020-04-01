#include <Parallel/MPI.h>
#include "GhostTimeCluster.h"

namespace seissol::time_stepping {
void GhostTimeCluster::sendCopyLayer(){
  SCOREP_USER_REGION( "sendCopyLayer", SCOREP_USER_REGION_TYPE_FUNCTION )

  for (unsigned int region = 0; region < meshStructure->numberOfRegions; ++region) {
    if (meshStructure->neighboringClusters[region][1] == static_cast<int>(globalClusterId)) {
      std::cout << seissol::MPI::mpi.rank() << " cluster " << globalClusterId << " region: "
        << region << " out of "
      << meshStructure->numberOfRegions
      << " sending copy layer for tag " <<  timeData + meshStructure->sendIdentifiers[region] << std::endl;

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
}
void GhostTimeCluster::receiveGhostLayer(){
  SCOREP_USER_REGION( "receiveGhostLayer", SCOREP_USER_REGION_TYPE_FUNCTION )
  for (unsigned int region = 0; region < meshStructure->numberOfRegions; ++region) {
    if (meshStructure->neighboringClusters[region][1] == static_cast<int>(globalClusterId) ) {
      std::cout << seissol::MPI::mpi.rank() << " cluster " << globalClusterId << " region: " <<
        region << " recv copy layer for tag " <<  timeData + meshStructure->receiveIdentifiers[region] << std::endl;
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
  // Always check for receives/send for quicker MPI progress.
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
  sendCopyLayer();
  //std::cout << "MPI recv AdvancedPredicitonTimeMessage" << std::endl;
}
void GhostTimeCluster::handleAdvancedCorrectionTimeMessage(const NeighborCluster&) {
  assert(testForGhostLayerReceives());
  receiveGhostLayer();
  //std::cout << "AdvancedCorr " << ct.maxTimeStepSize << " " << msg.time << std::endl;
}

}
