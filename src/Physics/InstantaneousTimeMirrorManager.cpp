// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "InstantaneousTimeMirrorManager.h"

#include "Equations/Datastructures.h"
#include "Initializer/CellLocalMatrices.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Model/CommonDatastructures.h"
#include "Modules/Module.h"
#include "Modules/Modules.h"
#include "SeisSol.h"

#include <cmath>
#include <cstddef>
#include <utils/logger.h>
#include <vector>

namespace seissol::ITM {

void InstantaneousTimeMirrorManager::init(double velocityScalingFactor,
                                          double triggerTime,
                                          seissol::geometry::MeshReader* meshReader,
                                          LTS::Storage& ltsStorage,
                                          const initializer::ClusterLayout* clusterLayout) {
  isEnabled_ = true; // This is to sync just before and after the ITM. This does not toggle the ITM.
                     // Need this by default as true for it to work.
  this->velocityScalingFactor_ = velocityScalingFactor;
  this->triggerTime_ = triggerTime;
  this->meshReader_ = meshReader;
  this->ltsStorage_ = &ltsStorage;
  this->clusterLayout_ = clusterLayout; // An empty timestepping is added. Need to discuss what
                                        // exactly is to be sent here
  setSyncInterval(triggerTime);
  Modules::registerHook(*this, ModuleHook::SynchronizationPoint);
}

void InstantaneousTimeMirrorManager::syncPoint(double currentTime) {
  Module::syncPoint(currentTime);

  logInfo() << "InstantaneousTimeMirrorManager: Factor " << velocityScalingFactor_;
  if (!isEnabled_) {
    logInfo() << "InstantaneousTimeMirrorManager: Skipping syncing at " << currentTime
              << "as it is disabled";
    return;
  }

  logInfo() << "InstantaneousTimeMirrorManager Syncing at " << currentTime;

  logInfo() << "Scaling velocitites by factor of " << velocityScalingFactor_;
  updateVelocities();

  logInfo() << "Updating CellLocalMatrices";
  initializer::initializeCellLocalMatrices(
      *meshReader_, *ltsStorage_, *clusterLayout_, seissolInstance_.getSeisSolParameters().model);
  // An empty timestepping is added. Need to discuss what exactly is to be sent here

  logInfo() << "Updating TimeSteps by a factor of " << 1 / velocityScalingFactor_;
  updateTimeSteps();

  logInfo() << "Finished flipping.";
  isEnabled_ = false;
}

template <typename MaterialType>
void InstantaneousTimeMirrorManager::updateVelocitiesForMaterialType() {
  logError() << "ITM material update is not implemented for this material type.";
}

template <>
void InstantaneousTimeMirrorManager::updateVelocitiesForMaterialType<seissol::model::ElasticMaterial>() {
  auto itmParameters = seissolInstance_.getSeisSolParameters().model.itmParameters;
  auto reflectionType = itmParameters.itmReflectionType;

  const auto updateMaterial = [&](model::Material& material) {
    const auto rho = material.getDensity();
    const auto lambda = material.getLambdaBar();
    const auto mu = material.getMuBar();

    if (reflectionType == seissol::initializer::parameters::ReflectionType::BothWaves) {
      material.setLameParameters(mu * velocityScalingFactor_ * velocityScalingFactor_,
                                 lambda * velocityScalingFactor_ * velocityScalingFactor_);
    } else if (reflectionType ==
               seissol::initializer::parameters::ReflectionType::BothWavesVelocity) {
      material.setDensity(rho * velocityScalingFactor_);
      material.setLameParameters(mu * velocityScalingFactor_, lambda * velocityScalingFactor_);
    } else if (reflectionType == seissol::initializer::parameters::ReflectionType::Pwave) {
      material.setLameParameters(mu, lambda * velocityScalingFactor_ * velocityScalingFactor_);
    } else if (reflectionType == seissol::initializer::parameters::ReflectionType::Swave) {
      const auto newLambda = (lambda + 2.0*mu)/velocityScalingFactor_ - 2.0*velocityScalingFactor_*mu;
      if (newLambda < 0.0) {
        logError() << "New lambda is negative. This is not allowed. Please adjust your scaling factor.";
      }      material.setDensity(velocityScalingFactor_*rho);
      material.setLameParameters(velocityScalingFactor_ * mu, newLambda);
    } else {
      logError() << "Unknown reflection type; material cannot be updated.";
    }
  };

  for (auto& layer : ltsStorage_->leaves(Ghost)) {
    auto* materialData = layer.var<LTS::MaterialData>();

#pragma omp parallel for schedule(static)
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      auto& material = materialData[cell];
      updateMaterial(material);
    }
  }
}

template <>
void InstantaneousTimeMirrorManager::updateVelocitiesForMaterialType<seissol::model::AnisotropicMaterial>() {
  auto itmParameters = seissolInstance_.getSeisSolParameters().model.itmParameters;

  for (auto& layer : ltsStorage_->leaves(Ghost)) {
    auto* materialData = layer.var<LTS::MaterialData>();

#pragma omp parallel for schedule(static)
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      auto& material = materialData[cell];
      material.setDensity(material.getDensity()/(velocityScalingFactor_*velocityScalingFactor_));
    }
  }
}

template <typename MaterialType>
void InstantaneousTimeMirrorManager::updateTimeStepsForMaterialType() {
  logError() << "ITM time-step update is not implemented for this material type.";
}

template <>
void InstantaneousTimeMirrorManager::updateTimeStepsForMaterialType<seissol::model::ElasticMaterial>() {
  auto itmParameters = seissolInstance_.getSeisSolParameters().model.itmParameters;
  auto reflectionType = itmParameters.itmReflectionType;

  if (reflectionType == seissol::initializer::parameters::ReflectionType::BothWaves ||
      reflectionType == seissol::initializer::parameters::ReflectionType::Pwave) {
    for (auto& cluster : clusters_) {
      cluster->setClusterTimes(cluster->getClusterTimes() / velocityScalingFactor_);
      auto* neighborClusters = cluster->getNeighborClusters();
      for (auto& neighborCluster : *neighborClusters) {
        neighborCluster.ct.setTimeStepSize(neighborCluster.ct.getTimeStepSize() /
                                           velocityScalingFactor_);
      }
    }
  }

  if (reflectionType == seissol::initializer::parameters::ReflectionType::Swave) {

    for (auto& cluster : clusters_) {
      cluster->setClusterTimes(cluster->getClusterTimes() * velocityScalingFactor_);
      auto* neighborClusters = cluster->getNeighborClusters();
      for (auto& neighborCluster : *neighborClusters) {
        neighborCluster.ct.setTimeStepSize(neighborCluster.ct.getTimeStepSize() *
                                           velocityScalingFactor_);
      }
    }
  }
}

template <>
void InstantaneousTimeMirrorManager::updateTimeStepsForMaterialType<seissol::model::AnisotropicMaterial>() {

  for (auto& cluster : clusters_) {
    cluster->setClusterTimes(cluster->getClusterTimes() / velocityScalingFactor_);
    auto* neighborClusters = cluster->getNeighborClusters();
    for (auto& neighborCluster : *neighborClusters) {
      neighborCluster.ct.setTimeStepSize(neighborCluster.ct.getTimeStepSize() /
                                         velocityScalingFactor_);
    }
  }
}

void InstantaneousTimeMirrorManager::updateVelocities() {
  updateVelocitiesForMaterialType<seissol::model::MaterialT>();
}

void InstantaneousTimeMirrorManager::updateTimeSteps() {
  updateTimeStepsForMaterialType<seissol::model::MaterialT>();
}

void InstantaneousTimeMirrorManager::setClusterVector(
    const std::vector<seissol::time_stepping::AbstractTimeCluster*>& clusters) {
  this->clusters_ = clusters;
}

void initializeTimeMirrorManagers(double scalingFactor,
                                  double triggerTime,
                                  seissol::geometry::MeshReader* meshReader,
                                  LTS::Storage& ltsStorage,
                                  InstantaneousTimeMirrorManager& increaseManager,
                                  InstantaneousTimeMirrorManager& decreaseManager,
                                  seissol::SeisSol& seissolInstance,
                                  const initializer::ClusterLayout* clusterLayout) {
  increaseManager.init(scalingFactor,
                       triggerTime,
                       meshReader,
                       ltsStorage,
                       clusterLayout); // An empty timestepping is added. Need to discuss what
                                       // exactly is to be sent here
  auto itmParameters = seissolInstance.getSeisSolParameters().model.itmParameters;
  const double eps = itmParameters.itmDuration;

  decreaseManager.init(1 / scalingFactor,
                       triggerTime + eps,
                       meshReader,
                       ltsStorage,
                       clusterLayout); // An empty timestepping is added. Need to discuss what
                                       // exactly is to be sent here
};
} // namespace seissol::ITM