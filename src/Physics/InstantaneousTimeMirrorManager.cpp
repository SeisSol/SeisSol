// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "InstantaneousTimeMirrorManager.h"
#include "Initializer/CellLocalMatrices.h"
#include "Modules/Modules.h"
#include "SeisSol.h"
#include <Initializer/Parameters/ModelParameters.h>
#include <Initializer/Typedefs.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/LTSTree.h>
#include <Memory/Tree/Layer.h>
#include <Memory/Tree/Lut.h>
#include <Modules/Module.h>
#include <cmath>
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

  logInfo() << "InstantaneousTimeMirrorManager: Factor " << velocityScalingFactor;
  if (!isEnabled) {
    logInfo() << "InstantaneousTimeMirrorManager: Skipping syncing at " << currentTime
              << "as it is disabled";
    return;
  }

  logInfo() << "InstantaneousTimeMirrorManager Syncing at " << currentTime;

  logInfo() << "Scaling velocitites by factor of " << velocityScalingFactor;
  updateVelocities();

  logInfo() << "Updating CellLocalMatrices";
  initializer::initializeCellLocalMatrices(*meshReader,
                                           ltsTree,
                                           lts,
                                           ltsLut,
                                           *timestepping,
                                           seissolInstance.getSeisSolParameters().model);
  // An empty timestepping is added. Need to discuss what exactly is to be sent here

  logInfo() << "Updating TimeSteps by a factor of " << 1 / velocityScalingFactor;
  updateTimeSteps();

  logInfo() << "Finished flipping.";
  isEnabled = false;
}

void InstantaneousTimeMirrorManager::updateVelocities() {
#ifdef USE_ANISOTROPIC
  logError() << "This feature has not been implemented for anisotropic yet";
#else
  auto itmParameters = seissolInstance.getSeisSolParameters().model.itmParameters;
  auto reflectionType = itmParameters.itmReflectionType;
  for (auto& layer : ltsTree->leaves(Ghost)) {
    CellMaterialData* materials = layer.var(lts->material);

    if (reflectionType == seissol::initializer::parameters::ReflectionType::BothWaves) {
      for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
        auto& material = materials[cell];
// Refocusing both waves
#ifndef USE_ACOUSTIC
        material.local.mu *= velocityScalingFactor * velocityScalingFactor;
#endif
        material.local.lambda *= velocityScalingFactor * velocityScalingFactor;
        for (int i = 0; i < 4; i++) {
#ifndef USE_ACOUSTIC
          material.neighbor[i].mu *= velocityScalingFactor * velocityScalingFactor;
#endif
          material.neighbor[i].lambda *= velocityScalingFactor * velocityScalingFactor;
        }
      }
    }

    if (reflectionType == seissol::initializer::parameters::ReflectionType::BothWavesVelocity) {
      for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
        auto& material = materials[cell];
        // Refocusing both waves with constant velocities
        material.local.lambda *= velocityScalingFactor;
#ifndef USE_ACOUSTIC
        material.local.mu *= velocityScalingFactor;
#endif
        material.local.rho *= velocityScalingFactor;
        for (int i = 0; i < 4; i++) {
          material.neighbor[i].lambda *= velocityScalingFactor;
#ifndef USE_ACOUSTIC
          material.neighbor[i].mu *= velocityScalingFactor;
#endif
          material.neighbor[i].rho *= velocityScalingFactor;
        }
      }
    }

    if (reflectionType == seissol::initializer::parameters::ReflectionType::Pwave) {
      for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
        auto& material = materials[cell];
        // Refocusing only P-waves
        material.local.lambda *= velocityScalingFactor * velocityScalingFactor;
        for (int i = 0; i < 4; i++) {
          material.neighbor[i].lambda *= velocityScalingFactor * velocityScalingFactor;
        }
      }
    }

    if (reflectionType == seissol::initializer::parameters::ReflectionType::Swave) {
      for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
        auto& material = materials[cell];
        // Refocusing only S-waves
        // material.local.lambda =
        //     -2.0 * velocityScalingFactor * material.local.mu +
        //     (material.local.lambda + 2.0 * material.local.mu) / velocityScalingFactor;
        // material.local.mu *= velocityScalingFactor;
        // material.local.rho *= velocityScalingFactor;
        material.local.rho = material.local.rho * material.local.lambda /
                             (material.local.lambda + 2.0 * material.local.getMuBar() -
                              2.0 * material.local.getMuBar() * velocityScalingFactor);
        material.local.lambda = material.local.lambda + 2.0 * material.local.getMuBar() -
                                2.0 * material.local.getMuBar() * velocityScalingFactor;
#ifndef USE_ACOUSTIC
        material.local.mu = velocityScalingFactor * material.local.mu;
#endif

        for (int i = 0; i < 4; i++) {
          // material.neighbor[i].lambda =
          //     -2.0 * velocityScalingFactor * material.neighbor[i].mu +
          //     (material.neighbor[i].lambda + 2.0 * material.neighbor[i].mu) /
          //     velocityScalingFactor;
          // material.neighbor[i].mu *= velocityScalingFactor;
          // material.neighbor[i].rho *= velocityScalingFactor;
          material.neighbor[i].rho =
              material.neighbor[i].rho * material.neighbor[i].lambda /
              (material.neighbor[i].lambda + 2.0 * material.neighbor[i].getMuBar() -
               2.0 * material.neighbor[i].getMuBar() * velocityScalingFactor);
          material.neighbor[i].lambda =
              material.neighbor[i].lambda + 2.0 * material.neighbor[i].getMuBar() -
              2.0 * material.neighbor[i].getMuBar() * velocityScalingFactor;
#ifndef USE_ACOUSTIC
          material.neighbor[i].mu = velocityScalingFactor * material.neighbor[i].mu;
#endif
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
      auto* neighborClusters = cluster->getNeighborClusters();
      for (auto& neighborCluster : *neighborClusters) {
        neighborCluster.ct.setTimeStepSize(neighborCluster.ct.getTimeStepSize() /
                                           velocityScalingFactor);
      }
    }

    for (auto& cluster : *ghostTimeClusters) {
      cluster->setClusterTimes(cluster->getClusterTimes() / velocityScalingFactor);
      auto* ghostNeighborClusters = cluster->getNeighborClusters();
      for (auto& neighborcluster : *ghostNeighborClusters) {
        neighborcluster.ct.setTimeStepSize(neighborcluster.ct.getTimeStepSize() /
                                           velocityScalingFactor);
      }
    }
  }

  if (reflectionType ==
      seissol::initializer::parameters::ReflectionType::Swave) { // refocusing only S-waves

    if (abs(timeStepScalingFactor - 1.0) < 1e-6) {
      timeStepScalingFactor = 0.5;
    }

    if (abs(timeStepScalingFactor - 0.5) < 1e-6) {
      timeStepScalingFactor = 2.0;
    }

    for (auto& cluster : *timeClusters) {
      cluster->setClusterTimes(cluster->getClusterTimes() * timeStepScalingFactor);
      auto* neighborClusters = cluster->getNeighborClusters();
      for (auto& neighborCluster : *neighborClusters) {
        neighborCluster.ct.setTimeStepSize(neighborCluster.ct.getTimeStepSize() *
                                           timeStepScalingFactor);
      }
    }

    for (auto& cluster : *ghostTimeClusters) {
      cluster->setClusterTimes(cluster->getClusterTimes() * timeStepScalingFactor);
      auto* ghostNeighborClusters = cluster->getNeighborClusters();
      for (auto& neighborcluster : *ghostNeighborClusters) {
        neighborcluster.ct.setTimeStepSize(neighborcluster.ct.getTimeStepSize() *
                                           timeStepScalingFactor);
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
