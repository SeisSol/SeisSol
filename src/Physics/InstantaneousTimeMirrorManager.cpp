#include "InstantaneousTimeMirrorManager.h"
#include "Modules/Modules.h"
#include "Initializer/CellLocalMatrices.h"
#include "SeisSol.h"

namespace seissol::ITM {

void InstantaneousTimeMirrorManager::init(double velocityScalingFactor,
                                          double triggerTime,
                                          MeshReader* meshReader,
                                          initializers::LTSTree* ltsTree,
                                          initializers::LTS* lts,
                                          initializers::Lut* ltsLut,
                                          TimeStepping* timestepping) {
  isEnabled = true; // This is to sync just before and after the ITM. This does not toggle the ITM.
                    // Need this by default as true for it to work.
  this->velocityScalingFactor = velocityScalingFactor;
  this->triggerTime = triggerTime;
  this->meshReader = meshReader;
  this->ltsTree = ltsTree;
  this->lts = lts;
  this->ltsLut = ltsLut;
  this->timestepping = timestepping; // An empty timestepping is added. Need to discuss what exactly
                                     // is to be sent here
  setSyncInterval(triggerTime);
  Modules::registerHook(*this, SYNCHRONIZATION_POINT);
}

void InstantaneousTimeMirrorManager::syncPoint(double currentTime) {
  Module::syncPoint(currentTime);

  const auto rank = MPI::mpi.rank();
  logInfo(rank) << "InstantaneousTimeMirrorManager: Factor " << velocityScalingFactor;
  if (!isEnabled) {
    logInfo(rank) << "InstantaneousTimeMirrorManager: Skipping syncing at " << currentTime
                  << "as it is disabled";
    return;
  }

  logInfo(rank) << "InstantaneousTimeMirrorManager Syncing at " << currentTime;

  logInfo(rank) << "Scaling velocitites by factor of " << velocityScalingFactor;
  updateVelocities();

  logInfo(rank) << "Updating CellLocalMatrices";
  initializers::initializeCellLocalMatrices(
      *meshReader, ltsTree, lts, ltsLut, *timestepping); // An empty timestepping is added. Need to
                                                         // discuss what exactly is to be sent here

   logInfo(rank) << "Updating TimeSteps by a factor of " << 1/velocityScalingFactor;
   updateTimeSteps();

  logInfo(rank) << "Finished flipping.";
  isEnabled = false;
}

void InstantaneousTimeMirrorManager::updateVelocities() {
  for (auto it = ltsTree->beginLeaf(initializers::LayerMask(Ghost)); it != ltsTree->endLeaf();
       ++it) {
    CellMaterialData* materials = it->var(lts->material);
    for (unsigned cell = 0; cell < it->getNumberOfCells(); ++cell) {
      auto& material = materials[cell];
      // material.local.mu = material.local.mu / velocityScalingFactor;
      // material.local.lambda = material.local.lambda / velocityScalingFactor;
      // material.local.rho = material.local.rho * velocityScalingFactor;
      // material.local.rho *= this->velocityScalingFactor * this->velocityScalingFactor;
      material.local.mu *= velocityScalingFactor*velocityScalingFactor;
      material.local.lambda *= velocityScalingFactor*velocityScalingFactor;
      for(int i=0; i<4; i++){
        material.neighbor[i].mu *= velocityScalingFactor*velocityScalingFactor;
        material.neighbor[i].lambda *= velocityScalingFactor*velocityScalingFactor;
      }
    }
  }
}

void InstantaneousTimeMirrorManager::updateTimeSteps() {

  for (auto& cluster : *timeClusters) {
    cluster->getClusterTimes().getTimeStepSize() =
        cluster->getClusterTimes().getTimeStepSize() / velocityScalingFactor;

    auto neighborClusters = cluster->getNeighborClusters();
    for (auto& neighborCluster : *neighborClusters) {
      neighborCluster.ct.getTimeStepSize() =
          neighborCluster.ct.getTimeStepSize() / velocityScalingFactor;
    }
  }

  for (auto& cluster : *ghostTimeClusters) {
    cluster->getClusterTimes().getTimeStepSize() =
        cluster->getClusterTimes().getTimeStepSize() / velocityScalingFactor;
    auto ghostNeighborClusters = cluster->getNeighborClusters();
    for (auto& neighborcluster : *ghostNeighborClusters) {
      neighborcluster.ct.getTimeStepSize() =
          neighborcluster.ct.getTimeStepSize() / velocityScalingFactor;
    }
  }
}

void InstantaneousTimeMirrorManager::setTimeClusterVector(
    std::vector<std::unique_ptr<seissol::time_stepping::TimeCluster>>* clusters) {
  timeClusters = clusters;
}

void InstantaneousTimeMirrorManager::setGhostClusterVector(
    std::vector<std::unique_ptr<seissol::time_stepping::GhostTimeCluster>>* clusters) {
  ghostTimeClusters = clusters;
}

void initializeTimeMirrorManagers(double scalingFactor,
                                  double triggerTime,
                                  MeshReader* meshReader,
                                  initializers::LTSTree* ltsTree,
                                  initializers::LTS* lts,
                                  initializers::Lut* ltsLut,
                                  InstantaneousTimeMirrorManager& increaseManager,
                                  InstantaneousTimeMirrorManager& decreaseManager,
                                  TimeStepping* timestepping) {
  increaseManager.init(scalingFactor,
                       triggerTime,
                       meshReader,
                       ltsTree,
                       lts,
                       ltsLut,
                       timestepping); // An empty timestepping is added. Need to discuss what
                                      // exactly is to be sent here
  auto& memoryManager = seissol::SeisSol::main.getMemoryManager();
  auto ITMParamaters = memoryManager.getITMParameters();
  const double eps = ITMParamaters->getITMTime(); // TODO(Lukas) Use CFL condition for this
  // const double eps = 1.0;
  decreaseManager.init(1 / scalingFactor,
                       triggerTime + eps,
                       meshReader,
                       ltsTree,
                       lts,
                       ltsLut,
                       timestepping); // An empty timestepping is added. Need to discuss what
                                      // exactly is to be sent here
};

} // namespace seissol::ITM