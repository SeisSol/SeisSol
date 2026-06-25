// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "InstantaneousTimeMirrorManager.h"

#include "Equations/Datastructures.h"
#include "Initializer/Model/CellLocalMatrices.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Modules/Module.h"
#include "Modules/Modules.h"
#include "SeisSol.h"

#include <cstddef>
#include <type_traits>
#include <utils/logger.h>
#include <vector>

namespace seissol::ITM {

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

void InstantaneousTimeMirrorManager::init(double velocityScalingFactor,
                                          double triggerTime,
                                          seissol::geometry::MeshReader* meshReader,
                                          LTS::Storage& ltsStorage,
                                          const initializer::ClusterLayout* clusterLayout) {
  constexpr bool IsElastic =
      std::is_same_v<seissol::model::MaterialT, seissol::model::ElasticMaterial>;
  constexpr bool IsAnisotropic =
      std::is_same_v<seissol::model::MaterialT, seissol::model::AnisotropicMaterial>;
  if constexpr (!IsElastic && !IsAnisotropic) {
    logError() << "ITM material update is not implemented for this material type.";
  }

  if constexpr (IsAnisotropic) {
    const auto reflectionType =
        seissolInstance_.getSeisSolParameters().model.itmParameters.itmReflectionType;
    if (!isAnisotropicReflectionTypeSupported(reflectionType)) {
      logError() << "Anisotropic materials cannot have Pwave, Swave, and BothWavesVelocity.";
    }
  }

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

template <typename MaterialType>
void InstantaneousTimeMirrorManager::updateVelocitiesForMaterialType() {
  constexpr bool IsElastic = std::is_same_v<MaterialType, seissol::model::ElasticMaterial>;
  constexpr bool IsAnisotropic = std::is_same_v<MaterialType, seissol::model::AnisotropicMaterial>;

  const auto itmParameters = seissolInstance_.getSeisSolParameters().model.itmParameters;
  const auto reflectionType = itmParameters.itmReflectionType;

  const auto updateMaterial = [&](MaterialType& material) {
    if constexpr (IsElastic) {
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
    }

    if constexpr (IsAnisotropic) {
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

template <typename MaterialType>
void InstantaneousTimeMirrorManager::updateTimeStepsForMaterialType() {
  constexpr bool IsElastic = std::is_same_v<MaterialType, seissol::model::ElasticMaterial>;
  constexpr bool IsAnisotropic = std::is_same_v<MaterialType, seissol::model::AnisotropicMaterial>;
  if constexpr (!IsElastic && !IsAnisotropic) {
    return;
  }

  const auto itmParameters = seissolInstance_.getSeisSolParameters().model.itmParameters;
  const auto reflectionType = itmParameters.itmReflectionType;

  double timeStepScaling = 1.0;
  if constexpr (IsElastic) {
    timeStepScaling = getElasticTimeStepScalingFactor(reflectionType, velocityScalingFactor_);
  } else {
    timeStepScaling = 1.0 / velocityScalingFactor_;
  }

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
