#include "InstantaneousTimeMirrorManager.h"
#include "Modules/Modules.h"
#include "Initializer/CellLocalMatrices.h"
#include "SeisSol.h"

namespace seissol::ITM {

void InstantaneousTimeMirrorManager::init(double velocityScalingFactor,
                                          double triggerTime,
                                          seissol::geometry::MeshReader* meshReader,
                                          initializers::LTSTree* ltsTree,
                                          initializers::LTS* lts,
                                          initializers::Lut* ltsLut,
                                          const TimeStepping* timestepping) {
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

  logInfo(rank) << "Updating TimeSteps by a factor of " << 1 / velocityScalingFactor;
  updateTimeSteps();

  logInfo(rank) << "Finished flipping.";
  isEnabled = false;
}

void InstantaneousTimeMirrorManager::updateVelocities() {
#ifdef USE_ANISOTROPIC
  logError() << "This feature has not been implemented for anisotropic yet";
#else
  auto itmParameters = seissol::SeisSol::main.getSeisSolParameters().itmParameters;
  auto reflectionType = itmParameters.reflectionType;
  for (auto it = ltsTree->beginLeaf(initializers::LayerMask(Ghost)); it != ltsTree->endLeaf();
       ++it) {
    CellMaterialData* materials = it->var(lts->material);

    if (reflectionType == seissol::initializer::parameters::ReflectionType::bothwaves) {
      for (unsigned cell = 0; cell < it->getNumberOfCells(); ++cell) {
        auto& material = materials[cell];
        // Refocusing both waves
        material.local.mu *= velocityScalingFactor * velocityScalingFactor;
        material.local.lambda *= velocityScalingFactor * velocityScalingFactor;
        for (int i = 0; i < 4; i++) {
          material.neighbor[i].mu *= velocityScalingFactor * velocityScalingFactor;
          material.neighbor[i].lambda *= velocityScalingFactor * velocityScalingFactor;
        }
      }
    }

    if (reflectionType == seissol::initializer::parameters::ReflectionType::bothwaves_velocity) {
      for (unsigned cell = 0; cell < it->getNumberOfCells(); ++cell) {
        auto& material = materials[cell];
        // Refocusing both waves with constant velocities
        material.local.lambda *= velocityScalingFactor;
        material.local.mu *= velocityScalingFactor;
        material.local.rho *= velocityScalingFactor;
        for (int i = 0; i < 4; i++) {
          material.neighbor[i].lambda *= velocityScalingFactor;
          material.neighbor[i].mu *= velocityScalingFactor;
          material.neighbor[i].rho *= velocityScalingFactor;
        }
      }
    }

    if (reflectionType == seissol::initializer::parameters::ReflectionType::pwave) {
      for (unsigned cell = 0; cell < it->getNumberOfCells(); ++cell) {
        auto& material = materials[cell];
        // Refocusing only P-waves
        material.local.lambda *= velocityScalingFactor * velocityScalingFactor;
        for (int i = 0; i < 4; i++) {
          material.neighbor[i].lambda *= velocityScalingFactor * velocityScalingFactor;
        }
      }
    }

    if (reflectionType == seissol::initializer::parameters::ReflectionType::swave) {
      for (unsigned cell = 0; cell < it->getNumberOfCells(); ++cell) {
        auto& material = materials[cell];
        // Refocusing only S-waves
        material.local.lambda =
            -2.0 * velocityScalingFactor * material.local.mu +
            (material.local.lambda + 2.0 * material.local.mu) / velocityScalingFactor;
        material.local.mu *= velocityScalingFactor;
        material.local.rho *= velocityScalingFactor;

        for (int i = 0; i < 4; i++) {
          material.neighbor[i].lambda =
              -2.0 * velocityScalingFactor * material.neighbor[i].mu +
              (material.neighbor[i].lambda + 2.0 * material.neighbor[i].mu) / velocityScalingFactor;
          material.neighbor[i].mu *= velocityScalingFactor;
          material.neighbor[i].rho *= velocityScalingFactor;
        }
      }
    }
  }
#endif
}

void InstantaneousTimeMirrorManager::updateTimeSteps() {

  auto itmParameters = seissol::SeisSol::main.getSeisSolParameters().itmParameters;
  auto reflectionType = itmParameters.reflectionType;

  if (reflectionType == seissol::initializer::parameters::ReflectionType::bothwaves ||
      reflectionType == seissol::initializer::parameters::ReflectionType::pwave)
  // refocusing both the waves. Default scenario. Works for both waves, only P-wave and constant
  // impedance case
  {
    for (auto& cluster : *timeClusters) {
      cluster->setClusterTimes(cluster->getClusterTimes() / velocityScalingFactor);
      auto neighborClusters = cluster->getNeighborClusters();
      for (auto& neighborCluster : *neighborClusters) {
        neighborCluster.ct.setTimeStepSize(neighborCluster.ct.getTimeStepSize() /
                                           velocityScalingFactor);
      }
    }

    for (auto& cluster : *ghostTimeClusters) {
      cluster->setClusterTimes(cluster->getClusterTimes() / velocityScalingFactor);
      auto ghostNeighborClusters = cluster->getNeighborClusters();
      for (auto& neighborcluster : *ghostNeighborClusters) {
        neighborcluster.ct.setTimeStepSize(neighborcluster.ct.getTimeStepSize() /
                                           velocityScalingFactor);
      }
    }
  }

  if (reflectionType ==
      seissol::initializer::parameters::ReflectionType::swave) { // refocusing only S-waves
    for (auto& cluster : *timeClusters) {
      //        cluster->getClusterTimes() = cluster->getClusterTimes() * velocityScalingFactor;
      cluster->setClusterTimes(cluster->getClusterTimes() * velocityScalingFactor);
      auto neighborClusters = cluster->getNeighborClusters();
      for (auto& neighborCluster : *neighborClusters) {
        neighborCluster.ct.setTimeStepSize(neighborCluster.ct.getTimeStepSize() *
                                           velocityScalingFactor);
      }
    }

    for (auto& cluster : *ghostTimeClusters) {
      //        cluster->getClusterTimes() = cluster->getClusterTimes() * velocityScalingFactor;
      cluster->setClusterTimes(cluster->getClusterTimes() * velocityScalingFactor);
      auto ghostNeighborClusters = cluster->getNeighborClusters();
      for (auto& neighborcluster : *ghostNeighborClusters) {
        neighborcluster.ct.setTimeStepSize(neighborcluster.ct.getTimeStepSize() *
                                           velocityScalingFactor);
      }
    }
  }
}

void InstantaneousTimeMirrorManager::setTimeClusterVector(
    std::vector<std::unique_ptr<seissol::time_stepping::TimeCluster>>* clusters) {
  timeClusters = clusters;
}

void InstantaneousTimeMirrorManager::setGhostClusterVector(
    std::vector<std::unique_ptr<seissol::time_stepping::AbstractGhostTimeCluster>>* clusters) {
  ghostTimeClusters = clusters;
}

void initializeTimeMirrorManagers(double scalingFactor,
                                  double triggerTime,
                                  seissol::geometry::MeshReader* meshReader,
                                  initializers::LTSTree* ltsTree,
                                  initializers::LTS* lts,
                                  initializers::Lut* ltsLut,
                                  InstantaneousTimeMirrorManager& increaseManager,
                                  InstantaneousTimeMirrorManager& decreaseManager,
                                  const TimeStepping* timestepping) {
  increaseManager.init(scalingFactor,
                       triggerTime,
                       meshReader,
                       ltsTree,
                       lts,
                       ltsLut,
                       timestepping); // An empty timestepping is added. Need to discuss what
                                      // exactly is to be sent here
  auto itmParameters = seissol::SeisSol::main.getSeisSolParameters().itmParameters;
  double eps = itmParameters.ITMTime;

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