// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "InstantaneousTimeMirrorManager.h"

#include "Initializer/CellLocalMatrices.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Model/CommonDatastructures.h"
#include "Modules/Module.h"
#include "Modules/Modules.h"
#include "SeisSol.h"

#include <cstddef>
#include <utils/logger.h>
#include <vector>

namespace seissol::physics {

bool isAnisotropicReflectionTypeSupported(
    seissol::initializer::parameters::ReflectionType reflectionType) {
  return reflectionType == seissol::initializer::parameters::ReflectionType::BothWaves;
}

double getSwaveScaledLambda(double lambda, double mu, double velocityScalingFactor) {
  return (lambda + 2.0 * mu) / velocityScalingFactor - 2.0 * velocityScalingFactor * mu;
}

double
    getElasticTimeStepScalingFactor(seissol::initializer::parameters::ReflectionType reflectionType,
                                    double velocityScalingFactor) {
  if (reflectionType == seissol::initializer::parameters::ReflectionType::BothWaves ||
      reflectionType == seissol::initializer::parameters::ReflectionType::Pwave) {
    return 1.0 / velocityScalingFactor;
  }
  if (reflectionType == seissol::initializer::parameters::ReflectionType::Swave) {
    return velocityScalingFactor;
  }
  return 1.0;
}

namespace {
void checkSupported(const model::Material& material,
                    initializer::parameters::ReflectionType reflectionType) {
  if (material.getMaterialType() != model::MaterialType::Elastic &&
      material.getMaterialType() != model::MaterialType::Anisotropic) {
    logError() << "ITM material update is not implemented for this material type.";
  }
  if (material.getMaterialType() == model::MaterialType::Anisotropic &&
      !isAnisotropicReflectionTypeSupported(reflectionType)) {
    logError() << "Anisotropic materials cannot have Pwave, Swave, and BothWavesVelocity.";
  }
}
} // namespace

InstantaneousTimeMirrorManager::InstantaneousTimeMirrorManager(seissol::SeisSol& seissolInstance)
    : seissolInstance_(seissolInstance) {}

InstantaneousTimeMirrorManager::~InstantaneousTimeMirrorManager() = default;

void InstantaneousTimeMirrorManager::init(double velocityScalingFactor,
                                          double triggerTime,
                                          seissol::geometry::MeshReader* meshReader,
                                          LTS::Storage& ltsStorage,
                                          const initializer::ClusterLayout* clusterLayout) {

  const auto itmParameters = seissolInstance_.getSeisSolParameters().model.itmParameters;
  const auto reflectionType = itmParameters.itmReflectionType;

  // check over all cells (cheap; though it can be reduced to at most one per layer)
  for (auto& layer : ltsStorage.leaves()) {
    for (std::size_t i = 0; i < layer.size(); ++i) {
      checkSupported(layer.cellRef(i).get<LTS::MaterialData>(), reflectionType);
    }
  }

  isEnabled_ = true; // This is to sync just before and after the ITM. This does not toggle the ITM.
                     // Need this by default as true for it to work.
  this->velocityScalingFactor_ = velocityScalingFactor;
  this->triggerTime_ = triggerTime;
  this->meshReader_ = meshReader;
  this->ltsStorage_ = &ltsStorage;
  this->clusterLayout_ = clusterLayout;
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

#ifdef ACL_DEVICE
  void* stream = device::DeviceInstance::getInstance().api->getDefaultStream();
  ltsStorage_->varSynchronizeTo<LTS::LocalIntegration>(
      seissol::initializer::AllocationPlace::Device, stream);
  ltsStorage_->varSynchronizeTo<LTS::NeighboringIntegration>(
      seissol::initializer::AllocationPlace::Device, stream);
  device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
#endif

  logInfo() << "Updating TimeSteps by a factor of " << 1 / velocityScalingFactor_;
  updateTimeSteps();

  logInfo() << "Finished flipping.";
  isEnabled_ = false;
}

void InstantaneousTimeMirrorManager::updateVelocities() {
  const auto itmParameters = seissolInstance_.getSeisSolParameters().model.itmParameters;
  const auto reflectionType = itmParameters.itmReflectionType;

  const auto updateMaterial = [&](model::Material& material) {
    if (material.getMaterialType() == model::MaterialType::Elastic) {
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
        const auto newLambda = getSwaveScaledLambda(lambda, mu, velocityScalingFactor_);
        if (newLambda < 0.0) {
          logError() << "New lambda is negative. This is not allowed. Please adjust your scaling "
                        "factor.";
        }
        material.setDensity(velocityScalingFactor_ * rho);
        material.setLameParameters(velocityScalingFactor_ * mu, newLambda);
      } else {
        logError() << "Unknown reflection type; material cannot be updated.";
      }
    } else if (material.getMaterialType() == model::MaterialType::Anisotropic) {
      // for anisotropic materials, you could scale down density
      // or scale up all the direction-dependent coefficients.
      // we scale density for code simplicity
      material.setDensity(material.getDensity() /
                          (velocityScalingFactor_ * velocityScalingFactor_));
    }
  };

  for (auto& layer : ltsStorage_->leaves(Ghost)) {

#pragma omp parallel for schedule(static)
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      auto& material = layer.cellRef(cell).get<LTS::MaterialData>();
      updateMaterial(material);
    }
  }
}

void InstantaneousTimeMirrorManager::updateTimeSteps() {
  const auto itmParameters = seissolInstance_.getSeisSolParameters().model.itmParameters;
  const auto reflectionType = itmParameters.itmReflectionType;

  const double timeStepScaling =
      getElasticTimeStepScalingFactor(reflectionType, velocityScalingFactor_);

  if (timeStepScaling != 1.0) {
    scaleClusterTimes(timeStepScaling);
  }
}

void InstantaneousTimeMirrorManager::scaleClusterTimes(double scalingFactor) {
  for (auto& cluster : clusters_) {
    cluster->setClusterTimes(cluster->getClusterTimes() * scalingFactor);
    auto* neighborClusters = cluster->getNeighborClusters();
    for (auto& neighborCluster : *neighborClusters) {
      neighborCluster.ct.setTimeStepSize(neighborCluster.ct.getTimeStepSize() * scalingFactor);
    }
  }
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
  increaseManager.init(scalingFactor, triggerTime, meshReader, ltsStorage, clusterLayout);
  auto itmParameters = seissolInstance.getSeisSolParameters().model.itmParameters;
  const double eps = itmParameters.itmDuration;

  decreaseManager.init(1 / scalingFactor, triggerTime + eps, meshReader, ltsStorage, clusterLayout);
};
} // namespace seissol::physics
