#include "InstantaneousTimeMirrorManager.h"
#include "Initializer/CellLocalMatrices.h"
#include "Modules/Modules.h"
#include "SeisSol.h"
#include <Initializer/LTS.h>
#include <Initializer/Parameters/ModelParameters.h>
#include <Initializer/tree/LTSTree.hpp>
#include <Initializer/tree/Layer.hpp>
#include <Initializer/tree/Lut.hpp>
#include <Initializer/typedefs.hpp>
#include <Modules/Module.h>
#include <memory>
#include <utils/logger.h>
#include <vector>

namespace seissol::ITM {

void InstantaneousTimeMirrorManager::init(double velocityScalingFactor,
                                          double triggerTime,
                                          seissol::geometry::MeshReader* meshReader,
                                          initializer::LTSTree* ltsTree,
                                          initializer::LTS* lts,
                                          initializer::Lut* ltsLut,
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
  Modules::registerHook(*this, ModuleHook::SynchronizationPoint);
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
  initializer::initializeCellLocalMatrices(
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
  auto itmParameters = seissolInstance.getSeisSolParameters().model.itmParameters;
  auto reflectionType = itmParameters.itmReflectionType;
  for (auto it = ltsTree->beginLeaf(initializer::LayerMask(Ghost)); it != ltsTree->endLeaf();
       ++it) {
    CellMaterialData* materials = it->var(lts->material);

    if (reflectionType == seissol::initializer::parameters::ReflectionType::BothWaves) {
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

    if (reflectionType == seissol::initializer::parameters::ReflectionType::BothWavesVelocity) {
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

    if (reflectionType == seissol::initializer::parameters::ReflectionType::Pwave) {
      for (unsigned cell = 0; cell < it->getNumberOfCells(); ++cell) {
        auto& material = materials[cell];
        // Refocusing only P-waves
        material.local.lambda *= velocityScalingFactor * velocityScalingFactor;
        for (int i = 0; i < 4; i++) {
          material.neighbor[i].lambda *= velocityScalingFactor * velocityScalingFactor;
        }
      }
    }

    if (reflectionType == seissol::initializer::parameters::ReflectionType::Swave) {
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

  auto itmParameters = seissolInstance.getSeisSolParameters().model.itmParameters;
  auto reflectionType = itmParameters.itmReflectionType;

  if (reflectionType == seissol::initializer::parameters::ReflectionType::BothWaves ||
      reflectionType == seissol::initializer::parameters::ReflectionType::Pwave)
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
      seissol::initializer::parameters::ReflectionType::Swave) { // refocusing only S-waves
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
                                  initializer::LTSTree* ltsTree,
                                  initializer::LTS* lts,
                                  initializer::Lut* ltsLut,
                                  InstantaneousTimeMirrorManager& increaseManager,
                                  InstantaneousTimeMirrorManager& decreaseManager,
                                  seissol::SeisSol& seissolInstance,
                                  const TimeStepping* timestepping) {
  increaseManager.init(scalingFactor,
                       triggerTime,
                       meshReader,
                       ltsTree,
                       lts,
                       ltsLut,
                       timestepping); // An empty timestepping is added. Need to discuss what
                                      // exactly is to be sent here
  auto itmParameters = seissolInstance.getSeisSolParameters().model.itmParameters;
  const double eps = itmParameters.itmDuration;

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
