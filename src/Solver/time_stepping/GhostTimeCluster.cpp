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
      std::cout << "test for receive succ" << std::endl;
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
      std::cout << "test for send succ" << std::endl;
      send = sendQueue.erase(send);
    } else {
      ++send;
    }
  }
  return sendQueue.empty();
}

bool GhostTimeCluster::act() {
  testForGhostLayerReceives();
  testForCopyLayerSends();
  bool yield = false;
  switch (state) {
    case ActorState::Corrected: {
      const auto minNeighbor = std::min_element(
          neighbors.begin(), neighbors.end(),
          [this] (NeighborCluster const& a, NeighborCluster const& b) {
            return a.ct.nextCorrectionTime(syncTime) < b.ct.nextCorrectionTime(syncTime);
          });
      const bool mayPredict = minNeighbor == neighbors.end() ||
          ct.predictionTime < minNeighbor->ct.nextCorrectionTime(syncTime) - timeTolerance;
      if (ct.correctionTime + timeTolerance >= syncTime) {
        std::cout << this << " " << globalClusterId << " SHOULD NOT HAPPEN " << 
          ct.correctionTime << " " << timeTolerance << " " << syncTime << std::endl;
        state = ActorState::Synced;
      } else if (testForGhostLayerReceives() && mayPredict) {
        std::cout << this << " " << "mpi pred dt_max=" << ct.maxTimeStepSize  << " dt=" << timeStepSize()
                  << " t_p=" << ct.predictionTime << " t_c=" << ct.correctionTime
                  << std::endl;
        ct.predictionTime += timeStepSize();
        for (auto& neighbor : neighbors) {
          std::cout << ct.maxTimeStepSize << " mpi send pred " << ct.predictionTime << " " << neighbor.ct.nextCorrectionTime(syncTime) <<
          " " << (ct.predictionTime >= neighbor.ct.nextCorrectionTime(syncTime)) << std::endl;
          if (ct.predictionTime >= neighbor.ct.nextCorrectionTime(syncTime) - timeTolerance)  {
            AdvancedPredictionTimeMessage message{};
            message.time = ct.predictionTime;
            neighbor.outbox->push(message);
          }
        }
        state = ActorState::Predicted;
      } else {
        yield = !processMessages();
      }

      break;
    }
    case ActorState::Predicted: {
      const auto minNeighbor = std::min_element(
          neighbors.begin(), neighbors.end(),
          [] (NeighborCluster const& a, NeighborCluster const& b) {
            return a.ct.predictionTime < b.ct.predictionTime;
          });
      const bool mayCorrect = minNeighbor == neighbors.end() ||
          ct.predictionTime <= minNeighbor->ct.predictionTime + timeTolerance;
      if (testForCopyLayerSends() && mayCorrect) {
        std::cout << "mpi corr dt_max=" << ct.maxTimeStepSize << " dt=" << timeStepSize()
                  << " t_p=" << ct.predictionTime << " t_c=" << ct.correctionTime
                  << std::endl;
        ct.correctionTime += timeStepSize();
        for (auto& neighbor : neighbors) {
          std::cout << ct.maxTimeStepSize << " mpi send corr " << ct.correctionTime << " " << neighbor.ct.predictionTime << std::endl;
          if (ct.correctionTime >= neighbor.ct.predictionTime - timeTolerance) {
            AdvancedCorrectionTimeMessage message{};
            message.time = ct.correctionTime;
            neighbor.outbox->push(message);
          }
        }
        state = ActorState::Corrected;
      } else {
        yield = !processMessages();
      }
      break;
    }
    case ActorState::Synced:
      if (ct.correctionTime + timeTolerance < syncTime) {
        std::cout << "trlalaalala " << ct.correctionTime << " " << syncTime << std::endl;
        receiveGhostLayer();
        state = ActorState::Corrected;
      } else {
        yield = true;
      }
      break;
  }
 // std::cout << "MPI " << ct.correctionTime << " at state " << static_cast<int>(state)
 // << " " << testForGhostLayerReceives() << std::endl;
  return yield;
}

bool GhostTimeCluster::processMessages() {
  bool processed = false;
  for (auto& neighbor : neighbors) {
    if (neighbor.inbox->hasMessages()) {
      processed = true;
      Message message = neighbor.inbox->pop();
      std::visit([&neighbor,this](auto&& msg) {
        using T = std::decay_t<decltype(msg)>;
        if constexpr (std::is_same_v<T, AdvancedPredictionTimeMessage>) {
          std::cout << "MPI recv AdvancedPredicitonTimeMessage" << std::endl;
          assert(msg.time > neighbor.ct.predictionTime);
          assert(testForCopyLayerSends());
          neighbor.ct.predictionTime = msg.time;
          sendCopyLayer();
          //std::cout << "AdvancedPred " << ct.maxTimeStepSize << " " << msg.time << std::endl;
        } else if constexpr (std::is_same_v<T, AdvancedCorrectionTimeMessage>) {
          assert(msg.time > neighbor.ct.correctionTime);
          neighbor.ct.correctionTime = msg.time;
          assert(testForGhostLayerReceives());
          receiveGhostLayer();
          //std::cout << "AdvancedCorr " << ct.maxTimeStepSize << " " << msg.time << std::endl;
        } else {
          static_assert(always_false<T>::value, "non-exhaustive visitor!");
        }
      }, message);
    }
  }
  return processed;
  }

}
